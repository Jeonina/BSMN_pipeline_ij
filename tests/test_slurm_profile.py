"""
SLURM cluster profile tests (Phase 3 — SPEC-BSMN-REFACTOR-001).

These tests verify:
1. workflow/profiles/slurm/config.yaml is valid and contains required keys
2. run.py exposes --cluster {none,slurm} and routes correctly
3. Helper: runtime minutes -> HH:MM:SS conversion

Snakemake version target: 8+ (executor-plugin-slurm).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import yaml

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROFILE_PATH = PROJECT_ROOT / "workflow" / "profiles" / "slurm" / "config.yaml"
RUN_PY = PROJECT_ROOT / "run.py"


# ─────────────────────────────────────────────────────────────────────────
# Profile YAML tests
# ─────────────────────────────────────────────────────────────────────────


def test_slurm_profile_yaml_is_valid() -> None:
    """Profile must exist and parse as a YAML mapping."""
    assert PROFILE_PATH.is_file(), f"SLURM profile not found: {PROFILE_PATH}"
    with PROFILE_PATH.open() as f:
        data = yaml.safe_load(f)
    assert isinstance(data, dict), "profile must be a YAML mapping"


def test_slurm_profile_required_keys() -> None:
    """Profile must declare jobs, restart-times, use-apptainer, and executor."""
    with PROFILE_PATH.open() as f:
        data = yaml.safe_load(f)

    required = ["jobs", "restart-times", "use-apptainer", "executor"]
    missing = [k for k in required if k not in data]
    assert not missing, f"profile missing required keys: {missing}"

    # Snakemake 8+ plugin name
    assert data["executor"] == "slurm", "executor must be 'slurm' (snakemake 8+ plugin)"
    assert isinstance(data["jobs"], int) and data["jobs"] >= 1
    assert isinstance(data["restart-times"], int) and data["restart-times"] >= 0
    assert data["use-apptainer"] is True


def test_slurm_profile_default_resources_map_correctly() -> None:
    """default-resources must reference snakemake's resource keys for SLURM mapping."""
    with PROFILE_PATH.open() as f:
        data = yaml.safe_load(f)

    assert "default-resources" in data, "profile must declare default-resources"
    defaults = data["default-resources"]
    # default-resources is a list of "key=value" strings in snakemake profile syntax,
    # or a mapping. Accept either form.
    if isinstance(defaults, list):
        joined = "\n".join(defaults)
    elif isinstance(defaults, dict):
        joined = "\n".join(f"{k}={v}" for k, v in defaults.items())
    else:
        pytest.fail(f"default-resources must be list or dict, got {type(defaults)}")

    # Must reference partition (the SLURM-specific key the plugin reads)
    assert "slurm_partition" in joined or "partition" in joined, (
        "default-resources must declare slurm_partition for the plugin"
    )


# ─────────────────────────────────────────────────────────────────────────
# run.py CLI tests
# ─────────────────────────────────────────────────────────────────────────


def test_run_py_accepts_cluster_slurm_flag() -> None:
    """`python run.py --help` output must advertise --cluster {none,slurm}."""
    result = subprocess.run(
        [sys.executable, str(RUN_PY), "--help"],
        capture_output=True,
        text=True,
        cwd=PROJECT_ROOT,
        check=False,
    )
    assert result.returncode == 0, f"--help failed: {result.stderr}"
    out = result.stdout
    assert "--cluster" in out, "run.py --help must list --cluster option"
    assert "slurm" in out and "none" in out, (
        f"--cluster must offer 'none' and 'slurm' choices; help text:\n{out}"
    )


def test_run_py_cluster_slurm_invokes_profile() -> None:
    """When --cluster slurm is given, the snakemake invocation uses --profile."""
    # Import lazily so the import failure (if any) surfaces as a test failure,
    # not a collection error.
    sys.path.insert(0, str(PROJECT_ROOT))
    import run as run_module  # noqa: E402

    captured: dict[str, list[str]] = {}

    def fake_run(cmd: list[str], *_args: object, **_kw: object) -> MagicMock:
        captured["cmd"] = list(cmd)
        m = MagicMock()
        m.returncode = 0
        return m

    # Patch input resolution to bypass sample TSV building.
    fake_rows = [{"sample_id": "s1", "readgroup": "rg", "fq1": "/x/r1.fq", "fq2": "/x/r2.fq"}]

    with (
        patch.object(run_module, "_resolve_inputs", return_value=fake_rows),
        patch.object(run_module, "write_tsv"),
        patch("subprocess.run", side_effect=fake_run),
        patch("shutil.which", return_value="/usr/bin/sbatch"),
        patch.object(sys, "argv", ["run.py", "/data/fq/", "--cluster", "slurm"]),
    ):
        with pytest.raises(SystemExit) as ei:
            run_module.main()
        assert ei.value.code == 0

    cmd = captured.get("cmd", [])
    assert "--profile" in cmd, f"snakemake cmd missing --profile: {cmd}"
    profile_idx = cmd.index("--profile")
    assert cmd[profile_idx + 1] == "workflow/profiles/slurm", (
        f"--profile value must be 'workflow/profiles/slurm', got {cmd[profile_idx + 1]}"
    )
    # When --cluster slurm is used, --cores should NOT be passed (profile handles concurrency)
    assert "--cores" not in cmd, f"--cores must not be passed with --cluster slurm: {cmd}"


def test_run_py_cluster_slurm_fails_without_sbatch() -> None:
    """If sbatch is not on PATH, run.py must exit non-zero with helpful message."""
    sys.path.insert(0, str(PROJECT_ROOT))
    import run as run_module  # noqa: E402

    fake_rows = [{"sample_id": "s1", "readgroup": "rg", "fq1": "/x/r1.fq", "fq2": "/x/r2.fq"}]

    captured_err: list[str] = []

    def fake_die(msg: str) -> None:
        captured_err.append(msg)
        raise SystemExit(1)

    with (
        patch.object(run_module, "_resolve_inputs", return_value=fake_rows),
        patch.object(run_module, "write_tsv"),
        patch.object(run_module, "_die", side_effect=fake_die),
        patch("shutil.which", return_value=None),
        patch.object(sys, "argv", ["run.py", "/data/fq/", "--cluster", "slurm"]),
    ):
        with pytest.raises(SystemExit) as ei:
            run_module.main()
        # Expect non-zero exit code AND error message about sbatch
        assert ei.value.code != 0
        assert captured_err, "expected _die() to be called with a SLURM/sbatch message"
        joined = " ".join(captured_err).lower()
        assert "sbatch" in joined or "slurm" in joined, (
            f"error message must mention sbatch or SLURM: {captured_err}"
        )


def test_run_py_no_cluster_flag_uses_local_cores() -> None:
    """Default behavior (no --cluster) must still pass --cores N to snakemake."""
    sys.path.insert(0, str(PROJECT_ROOT))
    import run as run_module  # noqa: E402

    captured: dict[str, list[str]] = {}

    def fake_run(cmd: list[str], *_a: object, **_k: object) -> MagicMock:
        captured["cmd"] = list(cmd)
        m = MagicMock()
        m.returncode = 0
        return m

    fake_rows = [{"sample_id": "s1", "readgroup": "rg", "fq1": "/x/r1.fq", "fq2": "/x/r2.fq"}]

    with (
        patch.object(run_module, "_resolve_inputs", return_value=fake_rows),
        patch.object(run_module, "write_tsv"),
        patch("subprocess.run", side_effect=fake_run),
        patch.object(sys, "argv", ["run.py", "/data/fq/", "--cores", "8"]),
    ):
        with pytest.raises(SystemExit) as ei:
            run_module.main()
        assert ei.value.code == 0

    cmd = captured.get("cmd", [])
    assert "--cores" in cmd, f"local mode must pass --cores: {cmd}"
    assert "--profile" not in cmd, f"local mode must NOT use --profile: {cmd}"
    cores_idx = cmd.index("--cores")
    assert cmd[cores_idx + 1] == "8", f"--cores value mismatch: {cmd[cores_idx + 1]}"


# ─────────────────────────────────────────────────────────────────────────
# Helper conversion
# ─────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "minutes,expected",
    [
        (30, "00:30:00"),
        (90, "01:30:00"),
        (1440, "24:00:00"),
        (0, "00:00:00"),
        (60, "01:00:00"),
        (1, "00:01:00"),
    ],
)
def test_runtime_minutes_to_hms_conversion(minutes: int, expected: str) -> None:
    """runtime_minutes_to_hms converts integer minutes to HH:MM:SS."""
    sys.path.insert(0, str(PROJECT_ROOT))
    import run as run_module  # noqa: E402

    assert hasattr(run_module, "runtime_minutes_to_hms"), (
        "run.py must define runtime_minutes_to_hms helper"
    )
    assert run_module.runtime_minutes_to_hms(minutes) == expected
