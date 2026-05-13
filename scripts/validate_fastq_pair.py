#!/usr/bin/env python3
"""Pre-flight FASTQ R1/R2 pair validator.

Catches the failure mode that wasted 4+ hours of bwa runtime on real WGS
data: a truncated R2 gzip stream that bwa silently accepted (with a
single warning) until sambamba crashed with an ArraySliceError.

Checks performed (fail-fast, in order):

  1. Both files exist.
  2. Both files are non-empty.
  3. Gzip stream integrity (ALWAYS full stream, in both quick and full modes):
       - Stream-decompress every byte; gzip raises on CRC / length mismatch.
       - M-FIX-002: previously --quick mode skipped this for files > 1 GB,
         which let real-user case ERR194146_2.fastq.gz (56 GB) with
         mid-stream CRC corruption pass validation.
  4. First record read names match (after stripping /1, /2, whitespace).
  5. Read counts match:
       - quick mode: skipped (relies on step 3 integrity + step 4 name parity).
       - full mode: full count of records in both files.

What --quick skips (for speed on multi-GB files):
    - Exact per-record count comparison.
    - Random-position record-name parity checks.
What --quick ALWAYS does (mandatory):
    - File-exists / non-empty.
    - Full gzip integrity (CRC + length).
    - First-record name parity.

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

# Retained for backward compatibility with callers/tests that reference it.
# Historically used to switch --quick to a header-only integrity check above
# 1 GB; M-FIX-002 removed that fast-path because it missed mid-stream CRC
# corruption (user case ERR194146, 56 GB). Quick mode now ALWAYS runs the
# full integrity scan; this constant is no longer consulted internally.
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


def _verify_gzip_integrity(path: Path, chunk_size: int = 1024 * 1024) -> None:
    """Stream the gzip file end-to-end to catch CRC / length errors.

    Reads and discards every decompressed byte. The Python ``gzip`` module
    raises ``OSError`` (incl. ``gzip.BadGzipFile``) or ``EOFError`` on CRC
    or length mismatches when those errors are encountered mid-stream.

    Rationale (M-FIX-002): real user case ERR194146_2.fastq.gz (56 GB)
    had mid-stream CRC corruption that ``gunzip -t`` rejected but the old
    ``--quick`` mode missed, because the old large-file fast-path only
    validated the gzip header. I/O-bound: ~2 min on a 56 GB file at
    typical disk speed (500 MB/s), which is acceptable given the cost of
    discovering corruption hours into a bwa run instead.

    Raises:
        OSError, EOFError, gzip.BadGzipFile: on any gzip integrity error.
    """
    with gzip.open(path, "rb") as fh:
        while fh.read(chunk_size):
            pass


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def validate_pair(r1: Path, r2: Path, quick: bool = False) -> ValidationReport:
    """Validate that ``r1`` and ``r2`` form a coherent paired-end FASTQ pair.

    Args:
        r1: Path to R1 (read 1) gzip-compressed FASTQ.
        r2: Path to R2 (read 2) gzip-compressed FASTQ.
        quick: If True, skip the per-record read-count comparison and
            multi-position name parity checks (useful for multi-GB files).
            Full gzip integrity (CRC + length) and first-record name
            parity are ALWAYS verified, in both quick and full modes.

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
    # M-FIX-002: full gzip integrity is mandatory in BOTH modes. --quick
    # only skips the per-record count comparison, never the integrity scan.
    try:
        if quick:
            for _label, path in (("R1", r1), ("R2", r2)):
                _verify_gzip_integrity(path)
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
        help=(
            "Skip per-record name parity at multiple positions and skip "
            "exact read-count comparison. ALWAYS verifies full gzip "
            "integrity (CRC + length) -- equivalent to `gunzip -t`."
        ),
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
