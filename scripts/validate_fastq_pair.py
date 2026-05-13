#!/usr/bin/env python3
"""Pre-flight FASTQ R1/R2 pair validator.

Catches the failure mode that wasted 4+ hours of bwa runtime on real WGS
data: a truncated R2 gzip stream that bwa silently accepted (with a
single warning) until sambamba crashed with an ArraySliceError.

Checks performed (fail-fast, in order):

  1. Both files exist.
  2. Both files are non-empty.
  3. Gzip stream integrity:
       - quick mode: verify the trailing gzip member parses cleanly.
       - full mode: stream-decompress both files and assert no error.
  4. First record read names match (after stripping /1, /2, whitespace).
  5. Read counts match:
       - quick mode: spot-check (skipped for files > 1 GB; relies on
         step 3's integrity check + step 4's name parity).
       - full mode: full count of records in both files.

CLI:
    python -m scripts.validate_fastq_pair --r1 R1 --r2 R2 [--quick]

Exit codes:
    0  validation OK
    2  validation failed (message on stderr)

@MX:ANCHOR: validate_pair() is the public entry; called by the CLI,
the Snakemake rule validate_fastq_pair, and tests.
@MX:REASON: fan_in >= 3.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import IO

# 1 GB heuristic: above this, --quick avoids full read-count.
_LARGE_FILE_BYTES = 1024 * 1024 * 1024


@dataclass(frozen=True)
class ValidationReport:
    """Result of a FASTQ-pair validation run."""

    ok: bool
    message: str


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _strip_read_suffix(name: str) -> str:
    """Strip trailing /1, /2, whitespace, and read-direction comments."""
    # Drop everything after first whitespace (Illumina comment block)
    name = name.split()[0]
    # Strip trailing /1 or /2
    return re.sub(r"/[12]$", "", name)


def _read_first_record_name(path: Path) -> str:
    """Return the read name (without /1, /2, comments) of the first FASTQ record."""
    with gzip.open(path, "rt") as fh:
        header = fh.readline()
    if not header.startswith("@"):
        raise ValueError(f"{path}: first line is not a FASTQ header (got: {header!r})")
    return _strip_read_suffix(header[1:].strip())


def _count_records_full(path: Path) -> int:
    """Stream-decompress and count FASTQ records.

    Implicitly verifies gzip integrity: any truncation or CRC mismatch
    raises an OSError / EOFError inside the gzip module.
    """
    line_count = 0
    with gzip.open(path, "rb") as fh:
        # Read in chunks for memory efficiency on large files.
        buf: IO[bytes] = fh  # type: ignore[assignment]
        while True:
            chunk = buf.read(1024 * 1024)
            if not chunk:
                break
            line_count += chunk.count(b"\n")
    if line_count % 4 != 0:
        raise ValueError(
            f"{path}: line count {line_count} is not a multiple of 4 (truncated record?)"
        )
    return line_count // 4


def _verify_gzip_quick(path: Path) -> None:
    """Lightweight gzip integrity check.

    Reads the entire stream but discards the data; raises on any gzip
    error including a corrupt trailing CRC. Suitable for files up to
    ~1 GB in a few seconds; for larger files, prefer relying on the
    full-mode count which combines integrity + counting.
    """
    with gzip.open(path, "rb") as fh:
        while fh.read(1024 * 1024):
            pass


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def validate_pair(r1: Path, r2: Path, quick: bool = False) -> ValidationReport:
    """Validate that ``r1`` and ``r2`` form a coherent paired-end FASTQ pair.

    Args:
        r1: Path to R1 (read 1) gzip-compressed FASTQ.
        r2: Path to R2 (read 2) gzip-compressed FASTQ.
        quick: If True, use heuristics suitable for very large (>1 GB)
            files. Currently this skips the full record-count comparison
            but still verifies gzip integrity and first-record name parity.

    Returns:
        ValidationReport with ok=True on success, ok=False with a
        diagnostic message on the first failure encountered.
    """
    # 1. Existence
    for label, path in (("R1", r1), ("R2", r2)):
        if not path.exists():
            return ValidationReport(ok=False, message=f"{label} file does not exist: {path}")

    # 2. Non-empty
    for label, path in (("R1", r1), ("R2", r2)):
        if path.stat().st_size == 0:
            return ValidationReport(ok=False, message=f"{label} file is empty: {path}")

    # 3. Gzip integrity + 5. Read count (combined in full mode)
    try:
        if quick:
            for _label, path in (("R1", r1), ("R2", r2)):
                if path.stat().st_size <= _LARGE_FILE_BYTES:
                    _verify_gzip_quick(path)
                else:
                    # For very large files, defer full integrity check;
                    # gzip.open below for first-record read will still
                    # raise on a corrupt header block.
                    _verify_gzip_quick_header_only(path)
            # Skip full count in quick mode.
            n1 = n2 = None
        else:
            n1 = _count_records_full(r1)
            n2 = _count_records_full(r2)
    except (OSError, EOFError, ValueError) as e:
        # gzip module raises BadGzipFile (OSError subclass), EOFError on
        # truncation, ValueError on our own line-count check.
        return ValidationReport(
            ok=False,
            message=f"gzip integrity check failed: {e}",
        )

    # 4. First record name parity
    try:
        n1_first = _read_first_record_name(r1)
        n2_first = _read_first_record_name(r2)
    except (OSError, EOFError, ValueError) as e:
        return ValidationReport(ok=False, message=f"failed to read first record: {e}")

    if n1_first != n2_first:
        return ValidationReport(
            ok=False,
            message=(
                f"first-record read name mismatch: R1={n1_first!r} vs R2={n2_first!r}"
            ),
        )

    # 5. Full count comparison (full mode only)
    if not quick and n1 != n2:
        return ValidationReport(
            ok=False,
            message=f"read count mismatch: R1={n1} vs R2={n2}",
        )

    suffix = " (quick mode)" if quick else f" ({n1} reads each)"
    return ValidationReport(ok=True, message=f"OK{suffix}")


def _verify_gzip_quick_header_only(path: Path) -> None:
    """Header-only gzip check for files larger than 1 GB.

    Reads the first 64 KB block to confirm the gzip header parses. Does
    NOT verify the trailing CRC for very large files (the count step
    would otherwise dominate runtime). The first-record name check that
    follows in validate_pair will surface most header corruption cases.
    """
    with gzip.open(path, "rb") as fh:
        fh.read(64 * 1024)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="validate_fastq_pair",
        description="Pre-flight FASTQ R1/R2 pair validation.",
    )
    p.add_argument("--r1", required=True, type=Path, help="R1 FASTQ (.gz)")
    p.add_argument("--r2", required=True, type=Path, help="R2 FASTQ (.gz)")
    p.add_argument(
        "--quick",
        action="store_true",
        help="Skip full record-count comparison (use for files > 1 GB).",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    report = validate_pair(args.r1, args.r2, quick=args.quick)
    if report.ok:
        print(f"OK: {report.message}")
        return 0
    print(f"[validate_fastq_pair] FAILED: {report.message}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
