#!/usr/bin/env python3
"""Validate config/samples.tsv against the REQ-SNK-007 schema.

The schema is binding (SPEC-BSMN-REFACTOR-001 §3.2 REQ-SNK-007 / AC-16):

  - tab-separated values
  - first row is a header containing **exactly** these four columns, in this
    order: ``sample_id``, ``readgroup``, ``fq1``, ``fq2``
  - at least one data row
  - ``fq1`` and ``fq2`` resolve to existing files (after ``~`` expansion
    and absolute-path resolution)

On any violation the validator exits with status 2 and writes a clear
diagnostic to stderr.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

EXPECTED_HEADER: tuple[str, ...] = ("sample_id", "readgroup", "fq1", "fq2")


@dataclass(frozen=True)
class SampleRecord:
    """One validated row of samples.tsv."""

    sample_id: str
    readgroup: str
    fq1: Path
    fq2: Path


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _die(message: str) -> None:
    """Print ``message`` to stderr and exit with status 2."""
    print(f"error: {message}", file=sys.stderr)
    raise SystemExit(2)


def _check_header(header: list[str]) -> None:
    if len(header) != len(EXPECTED_HEADER):
        _die(
            "samples.tsv header column count mismatch: "
            f"expected {len(EXPECTED_HEADER)} columns "
            f"({', '.join(EXPECTED_HEADER)}), got {len(header)}"
        )
    if tuple(header) != EXPECTED_HEADER:
        _die(
            "samples.tsv header must equal exactly: "
            f"{chr(9).join(EXPECTED_HEADER)} — got: {chr(9).join(header)}"
        )


def _resolve_fastq(raw: str, *, source: str) -> Path:
    if not raw:
        _die(f"samples.tsv {source}: empty path")
    path = Path(raw).expanduser()
    try:
        path = path.resolve(strict=False)
    except (OSError, RuntimeError) as exc:
        _die(f"samples.tsv {source}: cannot resolve path {raw!r}: {exc}")
    if not path.exists():
        _die(f"samples.tsv {source}: FASTQ file not found: {path}")
    return path


def _open_source(path: Path | str) -> tuple[TextIO, bool]:
    """Return (stream, should_close)."""
    if path == "-" or path == Path("-"):
        return sys.stdin, False
    p = Path(path)
    if not p.exists():
        _die(f"samples.tsv not found: {p}")
    return p.open("r", encoding="utf-8", newline=""), True


def _iter_data_rows(
    reader: Iterable[list[str]],
) -> list[list[str]]:
    return [row for row in reader if row and any(cell.strip() for cell in row)]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def validate(path: Path | str) -> list[SampleRecord]:
    """Validate ``path`` and return the parsed records.

    On failure, exits with status 2 via :func:`_die`.
    """
    stream, should_close = _open_source(path)
    try:
        reader = csv.reader(stream, delimiter="\t")
        rows = _iter_data_rows(reader)
    finally:
        if should_close:
            stream.close()

    if not rows:
        _die("samples.tsv is empty (no header, no data)")

    header = rows[0]
    _check_header(header)

    data_rows = rows[1:]
    if not data_rows:
        _die("samples.tsv has a header but no data rows")

    records: list[SampleRecord] = []
    for line_no, row in enumerate(data_rows, start=2):
        if len(row) != len(EXPECTED_HEADER):
            _die(
                f"samples.tsv line {line_no}: column count "
                f"{len(row)} does not match expected {len(EXPECTED_HEADER)} "
                f"(expected tab-separated)"
            )
        sample_id, readgroup, fq1_raw, fq2_raw = row
        if not sample_id.strip():
            _die(f"samples.tsv line {line_no}: sample_id is empty")
        if not readgroup.strip():
            _die(f"samples.tsv line {line_no}: readgroup is empty")
        fq1 = _resolve_fastq(fq1_raw, source=f"line {line_no} fq1")
        fq2 = _resolve_fastq(fq2_raw, source=f"line {line_no} fq2")
        records.append(
            SampleRecord(
                sample_id=sample_id.strip(),
                readgroup=readgroup.strip(),
                fq1=fq1,
                fq2=fq2,
            )
        )
    return records


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="validate_samples_tsv",
        description=(
            "Validate config/samples.tsv against the REQ-SNK-007 schema. "
            "Use '-' to read from stdin."
        ),
    )
    parser.add_argument(
        "path",
        help="Path to samples.tsv (or '-' to read from stdin)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    records = validate(args.path)
    print(
        f"OK: samples.tsv valid ({len(records)} sample(s))",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
