"""Tests for scripts.validate_fastq_pair (M-FIX-001 Bug 1).

Pre-flight FASTQ pair validation that runs BEFORE bwa_mem_sort to catch
R1/R2 mismatch, truncated gzip, and count mismatches in seconds rather
than hours.

@MX:ANCHOR: validate_pair is the only entry point used by the Snakemake
gating rule; signature changes break workflow/rules/mapping.smk.
@MX:REASON: fan_in >= 3 (CLI, Snakemake rule, future test cases).
"""

from __future__ import annotations

import gzip
import subprocess
import sys
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
