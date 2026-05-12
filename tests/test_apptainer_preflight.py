"""Tests for scripts.apptainer_preflight (REQ-SNK-008 / AC-17).

Required runtime floor: Apptainer >= 1.2.5 OR SingularityCE >= 3.11.0.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass

import pytest
from scripts import apptainer_preflight as a

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

@dataclass
class _FakeRun:
    """Stand-in for subprocess.CompletedProcess."""

    stdout: str = ""
    stderr: str = ""
    returncode: int = 0


def _make_runner(mapping: dict[str, _FakeRun | FileNotFoundError]):
    """Return a callable mimicking ``subprocess.run`` based on argv[0]."""

    def _run(cmd, *args, **kwargs):  # type: ignore[no-untyped-def]
        key = cmd[0] if isinstance(cmd, (list, tuple)) else cmd.split()[0]
        outcome = mapping.get(key)
        if outcome is None:
            raise FileNotFoundError(key)
        if isinstance(outcome, FileNotFoundError):
            raise outcome
        return subprocess.CompletedProcess(
            args=cmd,
            returncode=outcome.returncode,
            stdout=outcome.stdout,
            stderr=outcome.stderr,
        )

    return _run


# ---------------------------------------------------------------------------
# parse_version unit tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "raw,expected",
    [
        ("apptainer version 1.2.5", (1, 2, 5)),
        ("apptainer version 1.3.0-rc.1", (1, 3, 0)),
        ("singularity version 3.11.0", (3, 11, 0)),
        ("singularity-ce version 4.0.1", (4, 0, 1)),
        ("apptainer version 2.0.0", (2, 0, 0)),
    ],
)
def test_parse_version_happy(raw: str, expected: tuple[int, int, int]) -> None:
    assert a.parse_version(raw) == expected


@pytest.mark.parametrize("raw", ["", "not-a-version", "apptainer version foo"])
def test_parse_version_invalid_raises(raw: str) -> None:
    with pytest.raises(ValueError):
        a.parse_version(raw)


# ---------------------------------------------------------------------------
# preflight integration (subprocess.run is mocked)
# ---------------------------------------------------------------------------

def test_apptainer_meets_floor(monkeypatch: pytest.MonkeyPatch) -> None:
    runner = _make_runner(
        {"apptainer": _FakeRun(stdout="apptainer version 1.2.5\n")}
    )
    monkeypatch.setattr(a.subprocess, "run", runner)
    result = a.preflight()
    assert result.tool == "apptainer"
    assert result.version == (1, 2, 5)


def test_apptainer_above_floor(monkeypatch: pytest.MonkeyPatch) -> None:
    runner = _make_runner(
        {"apptainer": _FakeRun(stdout="apptainer version 1.4.0\n")}
    )
    monkeypatch.setattr(a.subprocess, "run", runner)
    result = a.preflight()
    assert result.version == (1, 4, 0)


def test_apptainer_below_floor_exits_2(monkeypatch: pytest.MonkeyPatch) -> None:
    runner = _make_runner(
        {"apptainer": _FakeRun(stdout="apptainer version 1.2.4\n")}
    )
    monkeypatch.setattr(a.subprocess, "run", runner)
    with pytest.raises(SystemExit) as exc:
        a.preflight()
    assert exc.value.code == 2


def test_singularity_fallback_meets_floor(monkeypatch: pytest.MonkeyPatch) -> None:
    runner = _make_runner(
        {
            "apptainer": FileNotFoundError("apptainer"),
            "singularity": _FakeRun(stdout="singularity version 3.11.0\n"),
        }
    )
    monkeypatch.setattr(a.subprocess, "run", runner)
    result = a.preflight()
    assert result.tool == "singularity"
    assert result.version == (3, 11, 0)


def test_singularity_below_floor_exits_2(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    runner = _make_runner(
        {
            "apptainer": FileNotFoundError("apptainer"),
            "singularity": _FakeRun(stdout="singularity version 3.10.5\n"),
        }
    )
    monkeypatch.setattr(a.subprocess, "run", runner)
    with pytest.raises(SystemExit) as exc:
        a.preflight()
    assert exc.value.code == 2


def test_neither_in_path_exits_2(monkeypatch: pytest.MonkeyPatch) -> None:
    runner = _make_runner(
        {
            "apptainer": FileNotFoundError("apptainer"),
            "singularity": FileNotFoundError("singularity"),
        }
    )
    monkeypatch.setattr(a.subprocess, "run", runner)
    with pytest.raises(SystemExit) as exc:
        a.preflight()
    assert exc.value.code == 2


def test_main_prints_ok(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    runner = _make_runner(
        {"apptainer": _FakeRun(stdout="apptainer version 1.2.5\n")}
    )
    monkeypatch.setattr(a.subprocess, "run", runner)
    rc = a.main([])
    captured = capsys.readouterr()
    assert rc == 0
    assert "OK" in captured.out
    assert "1.2.5" in captured.out
