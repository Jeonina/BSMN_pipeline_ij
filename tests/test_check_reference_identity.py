"""Tests for scripts.check_reference_identity (REQ-SNK-006 / AC-15).

Pins the reference selector to ``hg38_no_alt`` for all CI runs.
Aliases ``hg38_decoy``, ``hg38_v0``, and bare ``hg38`` MUST fail-fast.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from scripts import check_reference_identity as r

# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------

def test_env_var_correct(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("BSMN_REFERENCE", "hg38_no_alt")
    monkeypatch.delenv("BSMN_REFERENCE_FILE", raising=False)
    out = r.check()
    assert out == "hg38_no_alt"


def test_file_correct(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    marker = tmp_path / "ci_reference.txt"
    marker.write_text("hg38_no_alt\n", encoding="utf-8")
    monkeypatch.delenv("BSMN_REFERENCE", raising=False)
    monkeypatch.setenv("BSMN_REFERENCE_FILE", str(marker))
    out = r.check()
    assert out == "hg38_no_alt"


def test_env_var_wins_over_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    marker = tmp_path / "ci_reference.txt"
    marker.write_text("hg38_decoy\n", encoding="utf-8")
    monkeypatch.setenv("BSMN_REFERENCE", "hg38_no_alt")
    monkeypatch.setenv("BSMN_REFERENCE_FILE", str(marker))
    out = r.check()
    assert out == "hg38_no_alt"


# ---------------------------------------------------------------------------
# Failure paths
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "bad",
    ["hg38_decoy", "hg38_v0", "hg38", "GRCh38", ""],
)
def test_env_var_invalid_exits_2(
    monkeypatch: pytest.MonkeyPatch, bad: str
) -> None:
    monkeypatch.setenv("BSMN_REFERENCE", bad)
    monkeypatch.delenv("BSMN_REFERENCE_FILE", raising=False)
    with pytest.raises(SystemExit) as exc:
        r.check()
    assert exc.value.code == 2


def test_file_wrong_value_exits_2(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    marker = tmp_path / "ci_reference.txt"
    marker.write_text("hg38_decoy\n", encoding="utf-8")
    monkeypatch.delenv("BSMN_REFERENCE", raising=False)
    monkeypatch.setenv("BSMN_REFERENCE_FILE", str(marker))
    with pytest.raises(SystemExit) as exc:
        r.check()
    assert exc.value.code == 2


def test_no_source_configured_exits_2(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("BSMN_REFERENCE", raising=False)
    monkeypatch.delenv("BSMN_REFERENCE_FILE", raising=False)
    with pytest.raises(SystemExit) as exc:
        r.check()
    assert exc.value.code == 2


def test_file_missing_exits_2(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.delenv("BSMN_REFERENCE", raising=False)
    monkeypatch.setenv(
        "BSMN_REFERENCE_FILE", str(tmp_path / "does_not_exist.txt")
    )
    with pytest.raises(SystemExit) as exc:
        r.check()
    assert exc.value.code == 2
