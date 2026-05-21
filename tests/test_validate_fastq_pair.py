"""Tests for scripts.validate_fastq_pair (M-FIX-001 Bug 1, M-FIX-002).

Pre-flight FASTQ pair validation that runs BEFORE bwa_mem_sort to catch
R1/R2 mismatch, truncated gzip, and count mismatches in seconds rather
than hours.

M-FIX-002 (2026-05): user case ERR194146 R2 (56 GB) had mid-stream CRC
corruption that --quick mode previously missed. Quick mode now ALWAYS
performs full gzip integrity scan regardless of file size.

@MX:ANCHOR: validate_pair is the only entry point used by the Snakemake
gating rule; signature changes break workflow/rules/mapping.smk.
@MX:REASON: fan_in >= 3 (CLI, Snakemake rule, future test cases).
"""

from __future__ import annotations

import gzip
import subprocess
import sys
import time
from pathlib import Path

from scripts import validate_fastq_pair as v

# ---------------------------------------------------------------------------
# Helpers — synthetic FASTQ generators
# ---------------------------------------------------------------------------


def _synth_fastq(path: Path, n_reads: int, name_prefix: str = "READ") -> None:
    """Write a gzip-compressed FASTQ with n_reads synthetic records."""
    seq = "ACGT" * 25  # 100 bp
    qual = "I" * 100
    with gzip.open(path, "wt") as fh:
        for i in range(n_reads):
            fh.write(f"@{name_prefix}{i}\n{seq}\n+\n{qual}\n")


def _make_pair(tmp_path: Path, n: int = 100) -> tuple[Path, Path]:
    r1 = tmp_path / "sample_R1.fastq.gz"
    r2 = tmp_path / "sample_R2.fastq.gz"
    _synth_fastq(r1, n)
    _synth_fastq(r2, n)
    return r1, r2


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------


def test_validate_accepts_matched_synthetic_pair(tmp_path: Path) -> None:
    r1, r2 = _make_pair(tmp_path, n=100)
    report = v.validate_pair(r1, r2, quick=False)
    assert report.ok, f"Expected OK report, got: {report.message}"


# ---------------------------------------------------------------------------
# Failure paths
# ---------------------------------------------------------------------------


def test_validate_rejects_missing_file(tmp_path: Path) -> None:
    r1, _ = _make_pair(tmp_path, n=10)
    r2 = tmp_path / "does_not_exist.fastq.gz"
    report = v.validate_pair(r1, r2, quick=True)
    assert not report.ok
    assert "exist" in report.message.lower() or "not found" in report.message.lower()


def test_validate_rejects_empty_file(tmp_path: Path) -> None:
    r1, _ = _make_pair(tmp_path, n=10)
    r2 = tmp_path / "empty_R2.fastq.gz"
    r2.write_bytes(b"")
    report = v.validate_pair(r1, r2, quick=True)
    assert not report.ok
    assert "empty" in report.message.lower()


def test_validate_rejects_truncated_r2_gzip(tmp_path: Path) -> None:
    r1, r2 = _make_pair(tmp_path, n=100)
    # Corrupt R2 by truncating the trailer (last 16 bytes contain gzip CRC + ISIZE).
    raw = r2.read_bytes()
    r2.write_bytes(raw[: max(32, len(raw) - 1024)])
    report = v.validate_pair(r1, r2, quick=False)
    assert not report.ok
    assert "gzip" in report.message.lower()


def test_validate_rejects_count_mismatch(tmp_path: Path) -> None:
    r1 = tmp_path / "x_R1.fastq.gz"
    r2 = tmp_path / "x_R2.fastq.gz"
    _synth_fastq(r1, 100)
    _synth_fastq(r2, 50)
    report = v.validate_pair(r1, r2, quick=False)
    assert not report.ok
    assert "count" in report.message.lower()


def test_validate_rejects_first_record_name_mismatch(tmp_path: Path) -> None:
    r1 = tmp_path / "y_R1.fastq.gz"
    r2 = tmp_path / "y_R2.fastq.gz"
    _synth_fastq(r1, 10, name_prefix="A")
    _synth_fastq(r2, 10, name_prefix="B")
    report = v.validate_pair(r1, r2, quick=True)
    assert not report.ok
    msg = report.message.lower()
    assert "name" in msg or "mismatch" in msg


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def test_validate_cli_exits_with_correct_codes(tmp_path: Path) -> None:
    r1, r2 = _make_pair(tmp_path, n=20)

    # success → exit 0
    proc_ok = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.validate_fastq_pair",
            "--r1",
            str(r1),
            "--r2",
            str(r2),
            "--quick",
        ],
        capture_output=True,
        text=True,
        cwd=Path(__file__).resolve().parent.parent,
    )
    assert proc_ok.returncode == 0, proc_ok.stderr

    # missing file → exit 2
    proc_bad = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.validate_fastq_pair",
            "--r1",
            str(r1),
            "--r2",
            str(tmp_path / "missing.fastq.gz"),
            "--quick",
        ],
        capture_output=True,
        text=True,
        cwd=Path(__file__).resolve().parent.parent,
    )
    assert proc_bad.returncode == 2


def test_validate_quick_accepts_pair(tmp_path: Path) -> None:
    """--quick mode skips full count but still passes a valid pair."""
    r1, r2 = _make_pair(tmp_path, n=50)
    report = v.validate_pair(r1, r2, quick=True)
    assert report.ok


# ---------------------------------------------------------------------------
# M-FIX-002: --quick mode must catch mid-stream gzip corruption
# Regression: real user case ERR194146 R2 (56 GB) had mid-stream CRC error
# that --quick mode previously missed (only header was validated for >1 GB).
# ---------------------------------------------------------------------------


def _make_crc_corrupted_gzip(path: Path) -> None:
    """Create a gzip file with mid-stream CRC corruption.

    Writes ~100 synthetic FASTQ records (~10 KB compressed), then flips
    a byte after the gzip header to corrupt the deflate stream / CRC.
    """
    _synth_fastq(path, n_reads=100, name_prefix="READ")
    data = path.read_bytes()
    # Skip the first 18 bytes (gzip header) to avoid breaking the header.
    # Corrupt a byte near the middle of the stream.
    mid = len(data) // 2
    corrupt_pos = max(18, mid)
    corrupted = bytearray(data)
    corrupted[corrupt_pos] ^= 0xFF
    path.write_bytes(bytes(corrupted))


def test_quick_mode_rejects_mid_stream_crc_corruption(tmp_path: Path) -> None:
    """--quick mode MUST detect mid-stream CRC corruption (M-FIX-002).

    Real user case: ERR194146_2.fastq.gz (56 GB) had mid-stream CRC errors
    that gunzip -t rejected, but --quick mode previously passed because it
    only validated the gzip header on files > 1 GB.
    """
    r1 = tmp_path / "clean_R1.fastq.gz"
    r2 = tmp_path / "corrupt_R2.fastq.gz"
    _synth_fastq(r1, n_reads=100)
    _make_crc_corrupted_gzip(r2)

    report = v.validate_pair(r1, r2, quick=True)
    assert not report.ok, "Quick mode must reject mid-stream CRC corruption"
    msg = report.message.lower()
    assert "gzip" in msg or "crc" in msg or "integrity" in msg, (
        f"Expected gzip/CRC/integrity in message, got: {report.message}"
    )


def test_quick_mode_rejects_truncated_gzip(tmp_path: Path) -> None:
    """--quick mode MUST detect a gzip file truncated before its CRC trailer."""
    r1 = tmp_path / "clean_R1.fastq.gz"
    r2 = tmp_path / "truncated_R2.fastq.gz"
    _synth_fastq(r1, n_reads=100)
    _synth_fastq(r2, n_reads=500)  # bigger so truncation removes >1 KB of payload
    raw = r2.read_bytes()
    # Drop the last 1 KB, which removes the CRC + ISIZE trailer.
    r2.write_bytes(raw[: max(32, len(raw) - 1024)])

    report = v.validate_pair(r1, r2, quick=True)
    assert not report.ok, "Quick mode must reject truncated gzip"
    msg = report.message.lower()
    assert "gzip" in msg or "integrity" in msg or "eof" in msg, (
        f"Expected gzip/integrity in message, got: {report.message}"
    )


def test_quick_mode_cli_rejects_mid_stream_corruption(tmp_path: Path) -> None:
    """CLI in --quick mode must exit 2 on mid-stream gzip corruption."""
    r1 = tmp_path / "clean_R1.fastq.gz"
    r2 = tmp_path / "corrupt_R2.fastq.gz"
    _synth_fastq(r1, n_reads=100)
    _make_crc_corrupted_gzip(r2)

    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.validate_fastq_pair",
            "--r1",
            str(r1),
            "--r2",
            str(r2),
            "--quick",
        ],
        capture_output=True,
        text=True,
        cwd=Path(__file__).resolve().parent.parent,
    )
    assert proc.returncode == 2, f"Expected exit 2, got {proc.returncode}: {proc.stderr}"
    stderr_lower = proc.stderr.lower()
    assert "gzip" in stderr_lower or "crc" in stderr_lower or "integrity" in stderr_lower


def test_quick_mode_rejects_corruption_in_large_file_path(
    tmp_path: Path, monkeypatch: object
) -> None:
    """Regression for ERR194146: --quick mode on a "large" file path must
    still run a full integrity scan, NOT a header-only check.

    Reproduces the user case where R2 was 56 GB and triggered the
    large-file fast-path; we simulate by lowering the size threshold so a
    small test file is treated as "large".
    """
    import os as _os

    r1 = tmp_path / "clean_R1.fastq.gz"
    r2 = tmp_path / "corrupt_R2.fastq.gz"
    _synth_fastq(r1, n_reads=100)

    # Write random (incompressible) payload so the gzip output stays large
    # enough to hide corruption well past the first 64 KB read window.
    # ~1 MB random data -> ~1 MB compressed.
    with gzip.open(r2, "wb") as fh:
        fh.write(b"@READ0\n")
        fh.write(b"A" * 100 + b"\n+\n" + b"I" * 100 + b"\n")
        fh.write(_os.urandom(1024 * 1024))  # incompressible bulk

    raw = r2.read_bytes()
    assert len(raw) > 256 * 1024, (
        f"test file must exceed 64 KB header window (got {len(raw)} bytes)"
    )
    corrupted = bytearray(raw)
    # Corrupt at ~75% to ensure we are well past the 64 KB header window.
    corrupt_pos = (len(raw) * 3) // 4
    corrupted[corrupt_pos] ^= 0xFF
    r2.write_bytes(bytes(corrupted))

    # Force the >1 GB code path even though the file is ~hundreds of KB.
    monkeypatch.setattr(v, "_LARGE_FILE_BYTES", 100)  # type: ignore[attr-defined]

    report = v.validate_pair(r1, r2, quick=True)
    assert not report.ok, (
        "Quick mode on the large-file path must NOT skip mid-stream "
        "integrity (ERR194146 regression)"
    )
    msg = report.message.lower()
    assert "gzip" in msg or "crc" in msg or "integrity" in msg


def test_quick_mode_still_fast_on_clean_input(tmp_path: Path) -> None:
    """A ~10 MB clean gzip must validate in --quick mode quickly (< 5 s).

    Confirms the new full-stream integrity scan does not regress wall-clock
    cost on typical inputs. 10 MB of compressed FASTQ ~= ~200k reads.
    """
    r1 = tmp_path / "big_R1.fastq.gz"
    r2 = tmp_path / "big_R2.fastq.gz"
    # 200_000 reads * ~210 raw bytes per record ~= 42 MB raw -> ~10 MB gz.
    _synth_fastq(r1, n_reads=200_000)
    _synth_fastq(r2, n_reads=200_000)
    # Sanity: file is at least a few hundred KB so we are actually measuring
    # the scan (FASTQ with constant payload compresses very well).
    assert r2.stat().st_size > 500_000

    start = time.monotonic()
    report = v.validate_pair(r1, r2, quick=True)
    elapsed = time.monotonic() - start

    assert report.ok, f"Expected OK, got: {report.message}"
    assert elapsed < 5.0, f"Quick mode took {elapsed:.2f}s on 10 MB input (limit 5 s)"


# ---------------------------------------------------------------------------
# M-FIX-003: --quick mode must catch R1/R2 read-count mismatch
# Regression: real user case ERR194146 had R1 with N reads and R2 with N+k
# reads after fasterq-dump --split-files left singletons. The mismatch only
# surfaced 7.5 h into bwa_mem_sort as "paired reads have different names".
# Quick mode previously checked only the first record's name parity and let
# any mid/end-stream divergence slip through.
# ---------------------------------------------------------------------------


def test_quick_mode_rejects_r1_r2_count_mismatch(tmp_path: Path) -> None:
    """--quick mode MUST detect R1/R2 read-count inequality (M-FIX-003).

    Synthetic R1 with 100 reads, R2 with 99 reads. First-record name parity
    passes (both READ0). Without the parity count check, --quick mode would
    pass this pair and bwa would die hours later.
    """
    r1 = tmp_path / "uneven_R1.fastq.gz"
    r2 = tmp_path / "uneven_R2.fastq.gz"
    _synth_fastq(r1, n_reads=100)
    _synth_fastq(r2, n_reads=99)

    report = v.validate_pair(r1, r2, quick=True)
    assert not report.ok, "Quick mode must reject R1/R2 count mismatch"
    msg = report.message.lower()
    assert "count" in msg or "mismatch" in msg, (
        f"Expected count/mismatch in message, got: {report.message}"
    )


def test_quick_mode_rejects_invalid_fastq_structure(tmp_path: Path) -> None:
    """--quick mode MUST reject a gzip whose line count is not divisible by 4.

    A FASTQ record is exactly 4 lines (@header, seq, +, qual). A line count
    that is not a multiple of 4 implies a truncated or malformed record.
    """
    r1 = tmp_path / "clean_R1.fastq.gz"
    r2 = tmp_path / "malformed_R2.fastq.gz"
    # Both files have the SAME total line count so parity passes, but the
    # count is not a multiple of 4 -> structure check must fire.
    # 9 lines (one full 4-line record + a 5-line malformed second record).
    with gzip.open(r1, "wt") as fh:
        fh.write("@READ0\nACGT\n+\nIIII\n")
        fh.write("@READ1\nACGT\n+\nIIII\nEXTRA\n")  # 5 lines -> 9 total
    with gzip.open(r2, "wt") as fh:
        fh.write("@READ0\nACGT\n+\nIIII\n")
        fh.write("@READ1\nACGT\n+\nIIII\nEXTRA\n")  # 5 lines -> 9 total

    report = v.validate_pair(r1, r2, quick=True)
    assert not report.ok, "Quick mode must reject non-multiple-of-4 line count"
    msg = report.message.lower()
    assert "fastq" in msg or "structure" in msg or "multiple of 4" in msg, (
        f"Expected fastq/structure/multiple-of-4 in message, got: {report.message}"
    )


def test_quick_mode_accepts_matched_counts(tmp_path: Path) -> None:
    """--quick mode passes a pair with identical R1/R2 record counts (M-FIX-003)."""
    r1 = tmp_path / "matched_R1.fastq.gz"
    r2 = tmp_path / "matched_R2.fastq.gz"
    _synth_fastq(r1, n_reads=100)
    _synth_fastq(r2, n_reads=100)

    report = v.validate_pair(r1, r2, quick=True)
    assert report.ok, f"Expected OK on matched pair, got: {report.message}"
    assert "100" in report.message, f"Expected count in message, got: {report.message}"


def test_quick_mode_speed_budget(tmp_path: Path) -> None:
    """--quick mode on a ~1 MB pair completes well under 5 seconds (M-FIX-003).

    Proportional check: if 1 MB takes ~5 s, then a 60 GB pair extrapolates
    to ~85 hours, which would be unacceptable. Real budget on 1 MB should
    be sub-second; this guard exists to catch accidental quadratic regressions.
    """
    r1 = tmp_path / "speed_R1.fastq.gz"
    r2 = tmp_path / "speed_R2.fastq.gz"
    # ~20_000 reads * ~210 raw bytes ~= 4 MB raw, ~1 MB gz (with constant payload).
    _synth_fastq(r1, n_reads=20_000)
    _synth_fastq(r2, n_reads=20_000)
    assert r1.stat().st_size > 50_000, "test file must be large enough to time"

    start = time.monotonic()
    report = v.validate_pair(r1, r2, quick=True)
    elapsed = time.monotonic() - start

    assert report.ok, f"Expected OK, got: {report.message}"
    assert elapsed < 5.0, f"Quick mode took {elapsed:.2f}s on ~1 MB pair (budget 5 s)"
