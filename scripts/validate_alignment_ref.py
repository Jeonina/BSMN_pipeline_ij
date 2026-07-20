#!/usr/bin/env python3
"""Validate an external alignment's sequence dictionary against the reference.

Used by the ``ingest_alignment`` rule (workflow/rules/ingest.smk) to fail fast
when a pre-aligned BAM/CRAM was mapped to an incompatible reference. The check
compares the ``@SQ`` lines of the alignment header against the reference Picard
``.dict``:

  * A shared contig NAME with a different LENGTH is a hard failure (the
    references disagree on that contig and calling would be silently wrong).
  * A subset of reference contigs in the alignment is OK (e.g. a chr20-only
    BAM against a whole-genome reference).
  * Extra decoy/alt contigs present only in the alignment are OK.
  * Zero overlapping contig names means the references are unrelated → failure.

It also decides whether the output CRAM needs a read group injected: external
alignments may carry no ``@RG`` (e.g. the NHGRI novoalign HG002 BAM), or an
``@RG`` whose ``SM`` does not match the pipeline sample id. Mutect2's
``--tumor-sample {sample}`` matches against ``@RG SM`` tags, so a missing/
mismatched ``SM`` aborts calling with "samples cannot be empty". The
``--emit-rg-decision`` mode prints ``ok`` or ``inject`` for the ingest rule to
branch on.

This module is pure (string / file parsing) so it is unit-testable without
samtools or apptainer. The alignment header itself is produced on the cluster
via ``samtools view -H`` inside the samtools container and piped to ``--header``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def parse_sq_lines(header_text: str) -> dict[str, int]:
    """Parse ``@SQ`` lines from a SAM header into ``{contig_name: length}``.

    Lines without both ``SN:`` and ``LN:`` fields are ignored.
    """
    contigs: dict[str, int] = {}
    for line in header_text.splitlines():
        if not line.startswith("@SQ"):
            continue
        name: str | None = None
        length: int | None = None
        for field in line.split("\t")[1:]:
            if field.startswith("SN:"):
                name = field[3:]
            elif field.startswith("LN:"):
                try:
                    length = int(field[3:])
                except ValueError:
                    length = None
        if name is not None and length is not None:
            contigs[name] = length
    return contigs


def parse_dict_contigs(dict_path: Path) -> dict[str, int]:
    """Parse a Picard ``.dict`` file into ``{contig_name: length}``."""
    return parse_sq_lines(dict_path.read_text(encoding="utf-8"))


def check_compatibility(alignment: dict[str, int], reference: dict[str, int]) -> list[str]:
    """Return a list of human-readable error strings (empty == compatible).

    Compatibility rules are documented in the module docstring.
    """
    shared = set(alignment) & set(reference)
    if not shared:
        return [
            "no shared contigs between the alignment header and the reference "
            "dictionary; the BAM/CRAM appears to use an unrelated reference "
            f"(alignment contigs sample: {sorted(alignment)[:5]}, "
            f"reference contigs sample: {sorted(reference)[:5]})"
        ]

    errors: list[str] = []
    for name in sorted(shared):
        if alignment[name] != reference[name]:
            errors.append(
                f"contig '{name}' length mismatch: alignment LN={alignment[name]} "
                f"vs reference LN={reference[name]} (incompatible reference)"
            )
    return errors


def validate_against_dict(header_text: str, dict_path: Path) -> list[str]:
    """Parse both inputs and return compatibility errors (empty == OK)."""
    alignment = parse_sq_lines(header_text)
    if not alignment:
        return ["alignment header contains no @SQ lines (is it aligned?)"]
    reference = parse_dict_contigs(dict_path)
    if not reference:
        return [f"reference dictionary '{dict_path}' contains no @SQ lines"]
    return check_compatibility(alignment, reference)


def parse_rg_lines(header_text: str) -> list[dict[str, str]]:
    """Parse ``@RG`` lines from a SAM header into a list of tag dicts.

    Each returned dict maps two-letter tag codes to values, e.g.
    ``{"ID": "rg1", "SM": "HG002", "PL": "ILLUMINA"}``. Malformed fields
    without a ``KEY:VALUE`` shape are skipped.
    """
    read_groups: list[dict[str, str]] = []
    for line in header_text.splitlines():
        if not line.startswith("@RG"):
            continue
        tags: dict[str, str] = {}
        for field in line.split("\t")[1:]:
            key, sep, value = field.partition(":")
            if sep and key:
                tags[key] = value
        read_groups.append(tags)
    return read_groups


def decide_read_group(header_text: str, sample: str) -> str:
    """Decide whether the output CRAM needs a read group injected.

    Returns:
        ``"ok"``     — at least one ``@RG`` already has ``SM == sample``;
                       the existing read group(s) are usable for Mutect2's
                       ``--tumor-sample`` matching, so no rewrite is needed.
        ``"inject"`` — no ``@RG`` at all, OR an ``@RG`` exists but none has an
                       ``SM`` tag, OR no ``@RG`` has ``SM == sample``. The rule
                       must materialize a CRAM with a single read group whose
                       ``SM == sample`` so calling does not fail with
                       "samples cannot be empty".
    """
    read_groups = parse_rg_lines(header_text)
    for rg in read_groups:
        if rg.get("SM") == sample:
            return "ok"
    return "inject"


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="validate_alignment_ref",
        description=(
            "Validate an external alignment header. Two modes:\n"
            "  (default)            — fail (exit 1) if the @SQ lines are "
            "incompatible with the reference .dict.\n"
            "  --emit-rg-decision   — print 'ok' or 'inject' to stdout based on "
            "whether a read group with SM=={sample} already exists."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--header",
        required=True,
        help="Path to a file containing the SAM header (output of `samtools view -H`).",
    )
    parser.add_argument(
        "--dict",
        dest="dict_path",
        help="Path to the reference Picard .dict file (required unless --emit-rg-decision).",
    )
    parser.add_argument(
        "--emit-rg-decision",
        action="store_true",
        help="Print read-group decision ('ok'|'inject') to stdout and exit 0.",
    )
    parser.add_argument(
        "--sample",
        help="Sample id to match against @RG SM tags (required with --emit-rg-decision).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    header_text = Path(args.header).read_text(encoding="utf-8")

    # Mode 2: read-group decision (stdout = 'ok' | 'inject').
    if args.emit_rg_decision:
        if not args.sample:
            parser.error("--sample is required with --emit-rg-decision")
        print(decide_read_group(header_text, args.sample))
        return 0

    # Mode 1: reference-dictionary compatibility.
    if not args.dict_path:
        parser.error("--dict is required unless --emit-rg-decision is given")
    errors = validate_against_dict(header_text, Path(args.dict_path))
    if errors:
        print("error: alignment is incompatible with the reference:", file=sys.stderr)
        for err in errors:
            print(f"  - {err}", file=sys.stderr)
        return 1
    print("validate_alignment_ref: alignment is compatible with the reference.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
