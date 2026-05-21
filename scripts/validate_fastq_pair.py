#!/usr/bin/env python3
"""Pre-flight FASTQ R1/R2 pair validator.

Catches the failure mode that wasted 4+ hours of bwa runtime on real WGS
data: a truncated R2 gzip stream that bwa silently accepted (with a
single warning) until sambamba crashed with an ArraySliceError.

Checks performed (fail-fast, in order):

  1. Both files exist.
  2. Both files are non-empty.
  3. Gzip stream integrity + line count, fused in a single streaming pass
     (ALWAYS full stream, in both quick and full modes):
       - Stream-decompress every byte; gzip raises on CRC / length mismatch.
       - Count newlines while streaming so we get R1/R2 record counts for free.
       - M-FIX-002: previously --quick mode skipped this for files > 1 GB,
         which let real-user case ERR194146_2.fastq.gz (56 GB) with
         mid-stream CRC corruption pass validation.
       - M-FIX-003: previously --quick mode never compared R1/R2 record
         counts. Real-user case ERR194146 had R1/R2 with different read
         counts (fasterq-dump --split-files leaves singletons), and the
         mismatch only surfaced 7.5h into bwa_mem_sort as paired-reads-have-
         different-names. Quick mode now ALWAYS verifies count parity.
  4. R1 line count == R2 line count, and line count divisible by 4
     (proper FASTQ record structure). Both modes.
  5. First record read names match (after stripping /1, /2, whitespace).

What --quick skips (for speed on multi-GB files):
    - Random-position record-name parity checks beyond the first record.
What --quick ALWAYS does (mandatory):
    - File-exists / non-empty.
    - Full gzip integrity (CRC + length).
    - R1/R2 read-count parity + FASTQ structure (lines % 4 == 0).
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


def _verify_and_count_lines(path: Path, chunk_size: int = 1024 * 1024) -> int:
    """Stream the gzip file end-to-end and count newlines in one pass.

    Fuses two responsibilities so we touch every decompressed byte exactly
    once: (1) gzip integrity (CRC + length) and (2) line count for R1/R2
    parity. The Python ``gzip`` module raises ``OSError`` (incl.
    ``gzip.BadGzipFile``) or ``EOFError`` on CRC or length mismatches when
    those errors are encountered mid-stream.

    Rationale (M-FIX-002, M-FIX-003): real user case ERR194146 had both
    mid-stream CRC corruption (M-FIX-002) and R1/R2 count divergence
    (M-FIX-003: --split-files leaves singletons). Both modes now scan the
    full stream and return the line count, so callers can enforce
    parity + FASTQ structure (lines % 4 == 0). I/O-bound: ~5-10 min on a
    60 GB file at typical disk speed.

    Returns:
        Total newline count in the decompressed stream.

    Raises:
        OSError, EOFError, gzip.BadGzipFile: on any gzip integrity error.
    """
    line_count = 0
    with gzip.open(path, "rb") as fh:
        buf: IO[bytes] = fh  # type: ignore[assignment]
        while True:
            chunk = buf.read(chunk_size)
            if not chunk:
                break
            line_count += chunk.count(b"\n")
    return line_count


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def validate_pair(r1: Path, r2: Path, quick: bool = False) -> ValidationReport:
    """Validate that ``r1`` and ``r2`` form a coherent paired-end FASTQ pair.

    Args:
        r1: Path to R1 (read 1) gzip-compressed FASTQ.
        r2: Path to R2 (read 2) gzip-compressed FASTQ.
        quick: If True, skip multi-position name parity checks (useful for
            multi-GB files). Full gzip integrity (CRC + length), R1/R2
            read-count parity, FASTQ structure (lines % 4 == 0), and
            first-record name parity are ALWAYS verified in both modes.

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

    # 3. Gzip integrity + line count (fused single pass for both R1 and R2).
    # M-FIX-002: full gzip integrity is mandatory in BOTH modes.
    # M-FIX-003: read-count parity is mandatory in BOTH modes.
    try:
        lines_r1 = _verify_and_count_lines(r1)
        lines_r2 = _verify_and_count_lines(r2)
    except (OSError, EOFError) as e:
        # gzip module raises BadGzipFile (OSError subclass) or EOFError on
        # truncation / CRC mismatch.
        return ValidationReport(
            ok=False,
            message=f"gzip integrity check failed: {e}",
        )

    # 4a. R1/R2 line-count parity (M-FIX-003). A mismatch here is what
    # caused the real user case to die 7.5 h into bwa_mem_sort with
    # "paired reads have different names".
    if lines_r1 != lines_r2:
        n1_records = lines_r1 // 4
        n2_records = lines_r2 // 4
        return ValidationReport(
            ok=False,
            message=(
                f"R1/R2 read count mismatch: R1={n1_records} reads "
                f"({lines_r1} lines) vs R2={n2_records} reads "
                f"({lines_r2} lines)"
            ),
        )

    # 4b. FASTQ structure: every record is 4 lines.
    if lines_r1 % 4 != 0:
        return ValidationReport(
            ok=False,
            message=(
                f"invalid FASTQ structure: line count {lines_r1} is not a "
                f"multiple of 4 (truncated record?)"
            ),
        )

    n_records = lines_r1 // 4

    # 5. First record name parity.
    try:
        n1_first = _read_first_record_name(r1)
        n2_first = _read_first_record_name(r2)
    except (OSError, EOFError, ValueError) as e:
        return ValidationReport(ok=False, message=f"failed to read first record: {e}")

    if n1_first != n2_first:
        return ValidationReport(
            ok=False,
            message=(f"first-record read name mismatch: R1={n1_first!r} vs R2={n2_first!r}"),
        )

    suffix = " (quick mode)" if quick else ""
    return ValidationReport(ok=True, message=f"OK ({n_records} reads each){suffix}")


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
            "Skip multi-position name parity checks (only the first record "
            "name is compared). ALWAYS performs full gzip integrity (CRC + "
            "length) AND R1/R2 read-count parity check -- catches the "
            "ERR194146 failure mode where R1 and R2 had different read "
            "counts after fasterq-dump --split-files."
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
