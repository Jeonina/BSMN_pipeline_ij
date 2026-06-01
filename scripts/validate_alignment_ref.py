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


def check_compatibility(
    alignment: dict[str, int], reference: dict[str, int]
) -> list[str]:
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


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="validate_alignment_ref",
        description=(
            "Fail with a non-zero exit code if the alignment header @SQ lines "
            "are incompatible with the reference .dict."
        ),
    )
    parser.add_argument(
        "--header",
        required=True,
        help="Path to a file containing the SAM header (output of `samtools view -H`).",
    )
    parser.add_argument(
        "--dict",
        dest="dict_path",
        required=True,
        help="Path to the reference Picard .dict file.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    header_text = Path(args.header).read_text(encoding="utf-8")
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
