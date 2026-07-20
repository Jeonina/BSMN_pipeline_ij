"""CLI unit tests for the M2-consolidated filter scripts.

Per SPEC-BSMN-REFACTOR-001 §5 M2 + §3.3 REQ-DUP-002:

  - utils/somatic_vaf.2.py  → scripts/somatic_vaf_filter.py   (load_config stripped)
  - utils/strand_bias.2.py  → scripts/strand_bias_filter.py   (load_config stripped)
  - utils/repeat.2.py       → scripts/repeat_filter.py
  - utils/alt_bq_sum.py     → scripts/alt_bq_sum.py           (load_config stripped)

These tests verify:

  * the consolidated script file exists in scripts/,
  * its argparse signature is preserved against the winning *.2.py / .py base,
  * the legacy ``--reference`` / ``--conda-env`` flags are removed (load_config
    strip per REQ-DUP-002),
  * the module imports without raising — confirming no leftover
    ``from library.pileup import load_config`` reference.

Each test runs in well under 1s and uses no bio-tool subprocess.
Characterization (byte-equivalence) tests against M1b cluster goldens live
in ``test_char_*.py`` and skip until ``tests/data/golden/filtering/`` exists.
"""

from __future__ import annotations

import importlib
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def _has(mod: str) -> bool:
    return importlib.util.find_spec(mod) is not None


# Bio runtime dependencies (used by the production modules but not by
# every dev environment). Tests that require them skip when absent.
HAS_STATSMODELS = _has("statsmodels")
HAS_RPY2 = _has("rpy2")
HAS_SCIPY = _has("scipy")


def _run_help(script_rel: str) -> subprocess.CompletedProcess[str]:
    """Run ``python scripts/<script_rel> --help`` and return CompletedProcess."""
    return subprocess.run(
        [sys.executable, str(PROJECT_ROOT / "scripts" / script_rel), "--help"],
        capture_output=True,
        text=True,
        timeout=10,
    )


# =============================================================================
# scripts/somatic_vaf_filter.py
# =============================================================================


@pytest.mark.skipif(
    not HAS_STATSMODELS,
    reason="statsmodels not installed in this env (provided by Apptainer)",
)
class TestSomaticVafFilterCLI:
    SCRIPT = "somatic_vaf_filter.py"

    def test_script_file_exists(self) -> None:
        assert (PROJECT_ROOT / "scripts" / self.SCRIPT).is_file()

    def test_module_imports_without_load_config(self) -> None:
        mod = importlib.import_module("scripts.somatic_vaf_filter")
        # load_config MUST be stripped per REQ-DUP-002.
        assert not hasattr(mod, "load_config"), (
            "load_config should be stripped per REQ-DUP-002 — found leftover"
        )

    def test_help_succeeds(self) -> None:
        result = _run_help(self.SCRIPT)
        assert result.returncode == 0, result.stderr
        assert "--bam" in result.stdout
        assert "--min-MQ" in result.stdout or "--min-mq" in result.stdout.lower()
        assert "--nproc" in result.stdout

    def test_legacy_reference_and_conda_flags_removed(self) -> None:
        """REQ-DUP-002: ``--reference`` and ``--conda-env`` MUST be gone."""
        result = _run_help(self.SCRIPT)
        assert result.returncode == 0, result.stderr
        # The legacy load_config hook accepted --reference / --conda-env; both
        # are removed in the M2 consolidation.
        assert "--reference" not in result.stdout
        assert "--conda-env" not in result.stdout


# =============================================================================
# scripts/strand_bias_filter.py
# =============================================================================


@pytest.mark.skipif(
    not (HAS_RPY2 and HAS_SCIPY),
    reason="rpy2/scipy not installed in this env (provided by Apptainer)",
)
class TestStrandBiasFilterCLI:
    SCRIPT = "strand_bias_filter.py"

    def test_script_file_exists(self) -> None:
        assert (PROJECT_ROOT / "scripts" / self.SCRIPT).is_file()

    def test_module_imports_without_load_config(self) -> None:
        mod = importlib.import_module("scripts.strand_bias_filter")
        assert not hasattr(mod, "load_config")

    def test_help_succeeds(self) -> None:
        result = _run_help(self.SCRIPT)
        assert result.returncode == 0, result.stderr
        assert "--bam" in result.stdout
        assert "--nproc" in result.stdout

    def test_legacy_reference_and_conda_flags_removed(self) -> None:
        result = _run_help(self.SCRIPT)
        assert result.returncode == 0, result.stderr
        assert "--reference" not in result.stdout
        assert "--conda-env" not in result.stdout


# =============================================================================
# scripts/repeat_filter.py
# =============================================================================


class TestRepeatFilterCLI:
    SCRIPT = "repeat_filter.py"

    def test_script_file_exists(self) -> None:
        assert (PROJECT_ROOT / "scripts" / self.SCRIPT).is_file()

    def test_module_imports(self) -> None:
        # No load_config in the bare repeat.2.py to begin with — still verify
        # the module imports cleanly under the new location.
        importlib.import_module("scripts.repeat_filter")

    def test_help_succeeds(self) -> None:
        result = _run_help(self.SCRIPT)
        assert result.returncode == 0, result.stderr
        # Winner is repeat.2.py which uses ``--ref`` (not ``--reference``).
        assert "--ref" in result.stdout
        assert "--nproc" in result.stdout

    def test_no_library_pileup_reference(self) -> None:
        """Ensure no leftover ``from library.pileup`` import."""
        src = (PROJECT_ROOT / "scripts" / self.SCRIPT).read_text()
        assert "from library.pileup" not in src
        assert "from library." not in src


# =============================================================================
# scripts/alt_bq_sum.py
# =============================================================================


class TestAltBqSumCLI:
    """alt_bq_sum has no bio-runtime imports beyond bsmn_pipeline.pileup."""

    SCRIPT = "alt_bq_sum.py"

    def test_script_file_exists(self) -> None:
        assert (PROJECT_ROOT / "scripts" / self.SCRIPT).is_file()

    def test_module_imports_without_load_config(self) -> None:
        mod = importlib.import_module("scripts.alt_bq_sum")
        assert not hasattr(mod, "load_config"), (
            "load_config should be stripped per REQ-DUP-002 (DEF-002 fix)"
        )

    def test_help_succeeds(self) -> None:
        result = _run_help(self.SCRIPT)
        assert result.returncode == 0, result.stderr
        assert "--bam" in result.stdout
        assert "--min-BQ" in result.stdout or "--min-bq" in result.stdout.lower()

    def test_legacy_reference_and_conda_flags_removed(self) -> None:
        """DEF-002: ``--reference`` and ``--conda-env`` MUST be gone."""
        result = _run_help(self.SCRIPT)
        assert result.returncode == 0, result.stderr
        assert "--reference" not in result.stdout
        assert "--conda-env" not in result.stdout


# =============================================================================
# Legacy file deletion sanity (REQ-DUP-001)
# =============================================================================


class TestLegacyUtilsDeleted:
    """utils/ must be gone or, transitionally, contain none of these files."""

    LEGACY_PATHS = [
        "utils/somatic_vaf.py",
        "utils/somatic_vaf.2.py",
        "utils/strand_bias.py",
        "utils/strand_bias.2.py",
        "utils/repeat.py",
        "utils/repeat.2.py",
        "utils/PON_mask.2.py",
        "utils/germline_filter.py",
        "utils/alt_bq_sum.py",
        "utils/resubmit_hanging_jobs.sh",
    ]

    @pytest.mark.parametrize("legacy_path", LEGACY_PATHS)
    def test_legacy_file_removed(self, legacy_path: str) -> None:
        assert not (PROJECT_ROOT / legacy_path).exists(), (
            f"{legacy_path} should have been removed in M2 per REQ-DUP-001"
        )

    def test_utils_dir_removed(self) -> None:
        # AC-08: utils/ does not exist in the final tree.
        assert not (PROJECT_ROOT / "utils").is_dir(), (
            "utils/ directory should be removed in M2 per AC-08"
        )


# =============================================================================
# library/ collapse sanity (REQ-DUP-004 partial — job_queue retained for M7)
# =============================================================================


class TestLibraryCollapsed:
    """library/{config,misc,parser,pileup,__init__}.py are gone; job_queue.py stays."""

    GONE = [
        "library/__init__.py",
        "library/config.py",
        "library/misc.py",
        "library/parser.py",
        "library/pileup.py",
    ]

    @pytest.mark.parametrize("path", GONE)
    def test_legacy_library_module_removed(self, path: str) -> None:
        assert not (PROJECT_ROOT / path).exists(), f"{path} should be moved to bsmn_pipeline/ in M2"

    def test_job_queue_retained(self) -> None:
        # Held back to M7 per DEF-001 — jobs/*.py still depend on it.
        assert (PROJECT_ROOT / "library" / "job_queue.py").is_file()

    def test_bsmn_pipeline_modules_present(self) -> None:
        for name in ("config", "misc", "parser", "pileup"):
            assert (PROJECT_ROOT / "bsmn_pipeline" / f"{name}.py").is_file(), (
                f"bsmn_pipeline/{name}.py missing after M2 move"
            )
