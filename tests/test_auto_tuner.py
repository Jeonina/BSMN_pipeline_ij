"""
test_auto_tuner.py — Unit tests for the host-aware auto-tuner in
scripts/auto_params.py (derive_bwa_threads / derive_markdup /
derive_mutect2_threads / _tuning_params).

The derive_* functions are pure (cores, ram) -> knob, so most cases need no
mocking. _tuning_params queries the host, so its test patches os.cpu_count and
psutil to pin a synthetic host.

Evidence tiers under test:
    [DOC]       markdup thread ceiling 16; mutect2 default 4
    [HEURISTIC] bwa per-job width; >= 8-core spark cutoff
"""

from unittest.mock import MagicMock, patch

import pytest
from scripts.auto_params import (
    derive_bwa_threads,
    derive_markdup,
    derive_mutect2_threads,
    _tuning_params,
)

# Keys the workflow rules read out of resolved_params.yaml.
_REQUIRED_KEYS = {
    "bwa_threads",
    "sort_threads",
    "sort_memory",
    "markdup_engine",
    "markdup_threads",
    "mutect2_threads",
    "bqsr_memory_gb",
    "markdup_memory",
    "gatk_memory_gb",
}


class TestDeriveBwaThreads:
    """Per-job bwa width: aim ~12/job, bounded by RAM (~15 GB/job), share cores."""

    @pytest.mark.parametrize(
        ("cores", "ram_gb", "expected"),
        [
            (24, 62, 12),   # 62//15=4 ram-conc, 24//12=2 core-conc -> 2 jobs x 12
            (32, 128, 16),  # 128//15=8, 32//12=2 -> 2 jobs x 16 (matches storage32 intent)
            (180, 98, 30),  # 98//15=6 ram-bound, 15 core-conc -> 6 jobs x 30 (no idle cores)
            (8, 32, 8),     # small host -> one job on all cores
            (4, 16, 4),     # tiny host -> one job on all cores
            (64, 16, 64),   # RAM only fits 1 job -> one wide job, no oversubscription
        ],
    )
    def test_concurrency_aware(self, cores, ram_gb, expected):
        assert derive_bwa_threads(cores, ram_gb) == expected

    def test_never_below_one(self):
        assert derive_bwa_threads(0, 0) == 1
        assert derive_bwa_threads(1, 1) == 1

    def test_threads_times_concurrency_never_exceeds_cores(self):
        # A job must never be asked for more threads than the host has cores.
        for cores in (1, 2, 8, 16, 24, 32, 64, 128, 180):
            for ram in (8, 16, 32, 64, 128, 256):
                assert 1 <= derive_bwa_threads(cores, ram) <= cores


class TestDeriveMarkdup:
    """Auto-default is always Picard (spark is an explicit opt-in). The returned
    thread count (capped at 16 per DOC) applies only when spark is pinned."""

    @pytest.mark.parametrize(
        ("cores", "expected"),
        [
            (24, ("picard", 16)),   # engine always picard; spark-thread count capped at 16
            (16, ("picard", 16)),
            (8, ("picard", 8)),
            (7, ("picard", 7)),
            (4, ("picard", 4)),
            (200, ("picard", 16)),  # 16 ceiling holds far past 16
        ],
    )
    def test_engine_and_threads(self, cores, expected):
        assert derive_markdup(cores) == expected

    def test_engine_is_always_picard(self):
        for cores in (1, 4, 8, 16, 24, 168):
            assert derive_markdup(cores)[0] == "picard"

    def test_threads_never_exceed_doc_ceiling(self):
        for cores in (8, 16, 32, 64, 180):
            _, threads = derive_markdup(cores)
            assert threads <= 16


class TestDeriveMutect2Threads:
    """PairHMM threads: GATK default 4 (DOC), clamped down on tiny hosts."""

    @pytest.mark.parametrize(
        ("cores", "expected"),
        [(24, 4), (180, 4), (4, 4), (2, 2), (1, 1)],
    )
    def test_default_four_clamped_small(self, cores, expected):
        assert derive_mutect2_threads(cores) == expected


def _host(logical: int, physical: int, ram_gb: int):
    """Context managers pinning a synthetic host (cpu_count + psutil)."""
    mem = MagicMock()
    mem.total = ram_gb << 30
    mem.available = (ram_gb * 3 // 4) << 30
    return (
        patch("scripts.auto_params.os.cpu_count", return_value=logical),
        patch("scripts.auto_params.psutil.cpu_count", return_value=physical),
        patch("scripts.auto_params.psutil.virtual_memory", return_value=mem),
    )


class TestTuningParams:
    """_tuning_params emits every rule-facing key and honours the cores budget."""

    def test_emits_all_required_keys(self):
        cpu, cpuc, vm = _host(24, 12, 62)
        with cpu, cpuc, vm:
            params = _tuning_params()
        assert _REQUIRED_KEYS.issubset(params.keys())

    def test_cores_argument_overrides_detected_budget(self):
        # Machine has 180 logical cores, but the run was given --cores 24.
        cpu, cpuc, vm = _host(180, 90, 98)
        with cpu, cpuc, vm:
            params = _tuning_params(cores=24)
        assert params["tuned_for_cores"] == 24
        assert params["bwa_threads"] == derive_bwa_threads(24, 98)
        assert (params["markdup_engine"], params["markdup_threads"]) == derive_markdup(24)

    def test_default_uses_detected_logical_cpus(self):
        cpu, cpuc, vm = _host(24, 12, 62)
        with cpu, cpuc, vm:
            params = _tuning_params()
        assert params["tuned_for_cores"] == 24
        assert params["system_cpus"] == 24
        assert params["system_physical_cpus"] == 12

    def test_small_host_selects_picard(self):
        cpu, cpuc, vm = _host(4, 4, 16)
        with cpu, cpuc, vm:
            params = _tuning_params()
        assert params["markdup_engine"] == "picard"
