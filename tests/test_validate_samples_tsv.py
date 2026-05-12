"""Tests for scripts.validate_samples_tsv (REQ-SNK-007 / AC-16).

The samples.tsv schema is binding: exactly four tab-separated columns,
in this order, with a header row: sample_id, readgroup, fq1, fq2.
"""

from __future__ import annotations

import io
import subprocess
import sys
from pathlib import Path

import pytest
from scripts import validate_samples_tsv as v

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

VALID_HEADER = "sample_id\treadgroup\tfq1\tfq2"


def _write(tmp_path: Path, name: str, body: str) -> Path:
    p = tmp_path / name
    p.write_text(body, encoding="utf-8")
    return p


def _make_fastq(tmp_path: Path, name: str) -> Path:
    p = tmp_path / name
    p.write_bytes(b"\x1f\x8b\x08\x00")  # gzip magic prefix; content irrelevant
    return p


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------

def test_valid_tsv_returns_records(tmp_path: Path) -> None:
    fq1 = _make_fastq(tmp_path, "r1.fastq.gz")
    fq2 = _make_fastq(tmp_path, "r2.fastq.gz")
    tsv = _write(
        tmp_path,
        "samples.tsv",
        f"{VALID_HEADER}\nSAMPLE1\tRG1\t{fq1}\t{fq2}\n",
    )

    records = v.validate(tsv)

    assert len(records) == 1
    assert records[0].sample_id == "SAMPLE1"
    assert records[0].readgroup == "RG1"
    assert records[0].fq1 == fq1.resolve()
    assert records[0].fq2 == fq2.resolve()


def test_valid_tsv_multiple_rows(tmp_path: Path) -> None:
    fq1 = _make_fastq(tmp_path, "r1.fastq.gz")
    fq2 = _make_fastq(tmp_path, "r2.fastq.gz")
    tsv = _write(
        tmp_path,
        "samples.tsv",
        f"{VALID_HEADER}\nS1\tRG1\t{fq1}\t{fq2}\nS2\tRG2\t{fq1}\t{fq2}\n",
    )
    records = v.validate(tsv)
    assert len(records) == 2
    assert records[1].sample_id == "S2"


def test_tilde_path_is_expanded(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    fq1 = _make_fastq(tmp_path, "r1.fastq.gz")
    fq2 = _make_fastq(tmp_path, "r2.fastq.gz")
    # Point HOME at tmp_path so ~/r1.fastq.gz resolves to a real file.
    monkeypatch.setenv("HOME", str(tmp_path))
    tsv = _write(
        tmp_path,
        "samples.tsv",
        f"{VALID_HEADER}\nS\tRG\t~/r1.fastq.gz\t~/r2.fastq.gz\n",
    )
    records = v.validate(tsv)
    assert records[0].fq1.is_absolute()
    assert records[0].fq1 == fq1.resolve()
    assert records[0].fq2 == fq2.resolve()


def test_validate_from_stdin(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fq1 = _make_fastq(tmp_path, "r1.fastq.gz")
    fq2 = _make_fastq(tmp_path, "r2.fastq.gz")
    body = f"{VALID_HEADER}\nS\tRG\t{fq1}\t{fq2}\n"
    monkeypatch.setattr("sys.stdin", io.StringIO(body))
    records = v.validate("-")
    assert len(records) == 1
    assert records[0].sample_id == "S"


# ---------------------------------------------------------------------------
# Failure paths — each must raise SystemExit(2)
# ---------------------------------------------------------------------------

def test_missing_header_exits_2(tmp_path: Path) -> None:
    fq1 = _make_fastq(tmp_path, "r1.fastq.gz")
    fq2 = _make_fastq(tmp_path, "r2.fastq.gz")
    tsv = _write(tmp_path, "samples.tsv", f"S\tRG\t{fq1}\t{fq2}\n")
    with pytest.raises(SystemExit) as exc:
        v.validate(tsv)
    assert exc.value.code == 2


def test_wrong_column_order_exits_2(tmp_path: Path) -> None:
    fq1 = _make_fastq(tmp_path, "r1.fastq.gz")
    fq2 = _make_fastq(tmp_path, "r2.fastq.gz")
    # readgroup and sample_id swapped
    tsv = _write(
        tmp_path,
        "samples.tsv",
        f"readgroup\tsample_id\tfq1\tfq2\nRG\tS\t{fq1}\t{fq2}\n",
    )
    with pytest.raises(SystemExit) as exc:
        v.validate(tsv)
    assert exc.value.code == 2


def test_missing_column_exits_2(tmp_path: Path) -> None:
    fq1 = _make_fastq(tmp_path, "r1.fastq.gz")
    tsv = _write(
        tmp_path,
        "samples.tsv",
        f"sample_id\treadgroup\tfq1\nS\tRG\t{fq1}\n",
    )
    with pytest.raises(SystemExit) as exc:
        v.validate(tsv)
    assert exc.value.code == 2


def test_extra_column_exits_2(tmp_path: Path) -> None:
    fq1 = _make_fastq(tmp_path, "r1.fastq.gz")
    fq2 = _make_fastq(tmp_path, "r2.fastq.gz")
    tsv = _write(
        tmp_path,
        "samples.tsv",
        f"{VALID_HEADER}\textra\nS\tRG\t{fq1}\t{fq2}\tEXTRA\n",
    )
    with pytest.raises(SystemExit) as exc:
        v.validate(tsv)
    assert exc.value.code == 2


def test_empty_file_exits_2(tmp_path: Path) -> None:
    tsv = _write(tmp_path, "samples.tsv", "")
    with pytest.raises(SystemExit) as exc:
        v.validate(tsv)
    assert exc.value.code == 2


def test_comma_separated_exits_2(tmp_path: Path) -> None:
    fq1 = _make_fastq(tmp_path, "r1.fastq.gz")
    fq2 = _make_fastq(tmp_path, "r2.fastq.gz")
    tsv = _write(
        tmp_path,
        "samples.csv",
        f"sample_id,readgroup,fq1,fq2\nS,RG,{fq1},{fq2}\n",
    )
    with pytest.raises(SystemExit) as exc:
        v.validate(tsv)
    assert exc.value.code == 2


def test_fastq_not_found_exits_2(tmp_path: Path) -> None:
    tsv = _write(
        tmp_path,
        "samples.tsv",
        f"{VALID_HEADER}\nS\tRG\t{tmp_path}/missing_r1.fastq.gz\t{tmp_path}/missing_r2.fastq.gz\n",
    )
    with pytest.raises(SystemExit) as exc:
        v.validate(tsv)
    assert exc.value.code == 2


def test_no_data_rows_exits_2(tmp_path: Path) -> None:
    tsv = _write(tmp_path, "samples.tsv", f"{VALID_HEADER}\n")
    with pytest.raises(SystemExit) as exc:
        v.validate(tsv)
    assert exc.value.code == 2


# ---------------------------------------------------------------------------
# CLI behaviour
# ---------------------------------------------------------------------------

def test_cli_returns_zero_on_valid(tmp_path: Path) -> None:
    fq1 = _make_fastq(tmp_path, "r1.fastq.gz")
    fq2 = _make_fastq(tmp_path, "r2.fastq.gz")
    tsv = _write(
        tmp_path,
        "samples.tsv",
        f"{VALID_HEADER}\nS\tRG\t{fq1}\t{fq2}\n",
    )
    repo_root = Path(__file__).resolve().parent.parent
    result = subprocess.run(
        [sys.executable, "-m", "scripts.validate_samples_tsv", str(tsv)],
        cwd=repo_root,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr


def test_cli_returns_two_on_invalid(tmp_path: Path) -> None:
    tsv = _write(tmp_path, "samples.tsv", "")
    repo_root = Path(__file__).resolve().parent.parent
    result = subprocess.run(
        [sys.executable, "-m", "scripts.validate_samples_tsv", str(tsv)],
        cwd=repo_root,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
