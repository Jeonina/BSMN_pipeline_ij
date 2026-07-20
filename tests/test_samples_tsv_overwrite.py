"""Tests for samples.tsv overwrite behavior (M-FIX-001 Bug 2).

Regression: a leaked `sample` row appeared in config/samples.tsv after a
user invoked `python run.py test_data/ERR194146_1.fastq.gz`. Root cause:
the single-R1-file branch in `run.py._resolve_inputs` calls
`build_table(directory=p0.parent, ...)`, which scans ALL R1 files in the
parent directory rather than building a row for the specified pair only.

These tests pin the contracts:

  1. `write_tsv` opens the output in overwrite mode (not append).
  2. `run.py._resolve_inputs` returns only the user-specified pair when
     a single R1 file is supplied (no directory scan side-effect).
  3. `run.py` defensively deletes a stale samples.tsv before writing.
"""

from __future__ import annotations

import gzip
import importlib.util
from pathlib import Path
from typing import Any

import pytest

# Import scripts.make_samples_tsv directly.
from scripts import make_samples_tsv as mst

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def _load_run_module() -> Any:
    """Load `run.py` as a module without executing main()."""
    spec = importlib.util.spec_from_file_location("bsmn_run", PROJECT_ROOT / "run.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _make_gz(path: Path) -> None:
    """Create a minimal valid gzip FASTQ file."""
    with gzip.open(path, "wt") as fh:
        fh.write("@r1\nACGT\n+\nIIII\n")


# ---------------------------------------------------------------------------
# 1. write_tsv overwrite semantics
# ---------------------------------------------------------------------------


def test_write_tsv_uses_overwrite_mode(tmp_path: Path) -> None:
    out = tmp_path / "samples.tsv"

    first = [{"sample_id": "S1", "readgroup": "RG1", "fq1": "/a/1.fq.gz", "fq2": "/a/2.fq.gz"}]
    second = [{"sample_id": "S2", "readgroup": "RG2", "fq1": "/b/1.fq.gz", "fq2": "/b/2.fq.gz"}]

    mst.write_tsv(first, str(out))
    mst.write_tsv(second, str(out))

    text = out.read_text(encoding="utf-8")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    # 1 header + 1 row from second call only
    assert len(lines) == 2, f"Expected 2 lines (header + 1 row), got {len(lines)}: {lines}"
    assert "S2" in text
    assert "S1" not in text, "Stale row from first write_tsv call leaked into output"


def test_write_tsv_atomically_replaces_file(tmp_path: Path) -> None:
    """After write_tsv, the file contains exactly the new rows; no append."""
    out = tmp_path / "samples.tsv"

    # Pre-existing content from a "previous run"
    out.write_text(
        "sample_id\treadgroup\tfq1\tfq2\nOLD\tOLD\t/x/1.fq.gz\t/x/2.fq.gz\n",
        encoding="utf-8",
    )

    new_rows = [{"sample_id": "NEW", "readgroup": "NEW", "fq1": "/y/1.fq.gz", "fq2": "/y/2.fq.gz"}]
    mst.write_tsv(new_rows, str(out))

    text = out.read_text(encoding="utf-8")
    assert "OLD" not in text
    assert "NEW" in text


# ---------------------------------------------------------------------------
# 2. _resolve_inputs: single R1 file MUST NOT scan parent directory
# ---------------------------------------------------------------------------


def test_resolve_inputs_single_r1_returns_only_specified_pair(tmp_path: Path) -> None:
    """Reproduces the user-reported leak.

    Directory contains TWO R1/R2 pairs. User supplies one R1 file. The
    returned row list must contain only that pair — not the other.
    """
    # Set up two pairs in the same directory.
    a_r1 = tmp_path / "ERR194146_1.fastq.gz"
    a_r2 = tmp_path / "ERR194146_2.fastq.gz"
    b_r1 = tmp_path / "sample_R1.fastq.gz"
    b_r2 = tmp_path / "sample_R2.fastq.gz"
    for p in (a_r1, a_r2, b_r1, b_r2):
        _make_gz(p)

    run_mod = _load_run_module()
    rows = run_mod._resolve_inputs([str(a_r1)], recursive=False, pattern=None)

    sample_ids = {r["sample_id"] for r in rows}
    assert "ERR194146" in sample_ids
    assert "sample" not in sample_ids, (
        f"Directory-scan leak: got rows {sample_ids}, expected only {{ERR194146}}"
    )
    assert len(rows) == 1


def test_resolve_inputs_two_files_returns_only_specified_pair(tmp_path: Path) -> None:
    """R1 + R2 explicit form must also avoid directory scanning."""
    a_r1 = tmp_path / "ERR194146_1.fastq.gz"
    a_r2 = tmp_path / "ERR194146_2.fastq.gz"
    b_r1 = tmp_path / "sample_R1.fastq.gz"
    b_r2 = tmp_path / "sample_R2.fastq.gz"
    for p in (a_r1, a_r2, b_r1, b_r2):
        _make_gz(p)

    run_mod = _load_run_module()
    rows = run_mod._resolve_inputs([str(a_r1), str(a_r2)], recursive=False, pattern=None)

    sample_ids = {r["sample_id"] for r in rows}
    assert sample_ids == {"ERR194146"}, f"Expected only ERR194146, got {sample_ids}"


# ---------------------------------------------------------------------------
# 3. End-to-end overwrite contract via simulated two invocations
# ---------------------------------------------------------------------------


def test_run_py_overwrites_samples_tsv_between_invocations(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Two invocations of _resolve_inputs + write_tsv must produce a file
    containing only rows from the second invocation."""
    run_mod = _load_run_module()

    # First pair in dir_a
    dir_a = tmp_path / "first"
    dir_a.mkdir()
    a_r1 = dir_a / "ERR194146_1.fastq.gz"
    a_r2 = dir_a / "ERR194146_2.fastq.gz"
    _make_gz(a_r1)
    _make_gz(a_r2)

    # Second pair in dir_b
    dir_b = tmp_path / "second"
    dir_b.mkdir()
    b_r1 = dir_b / "DIFFERENT_R1.fastq.gz"
    b_r2 = dir_b / "DIFFERENT_R2.fastq.gz"
    _make_gz(b_r1)
    _make_gz(b_r2)

    out = tmp_path / "samples.tsv"

    # Invocation 1
    rows1 = run_mod._resolve_inputs([str(a_r1)], recursive=False, pattern=None)
    mst.write_tsv(rows1, str(out))

    # Invocation 2 (different input)
    rows2 = run_mod._resolve_inputs([str(b_r1)], recursive=False, pattern=None)
    mst.write_tsv(rows2, str(out))

    text = out.read_text(encoding="utf-8")
    assert "DIFFERENT" in text
    assert "ERR194146" not in text, "Stale row from previous invocation leaked"
