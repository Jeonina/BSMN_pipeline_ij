"""Unit tests for the pure planning/summary helpers in scripts/benchmark_bwa.py.
(The measurement path needs bwa + apptainer, so it is exercised on a host, not here.)"""

from scripts.benchmark_bwa import plan_runs, summarize


class TestPlanRuns:
    def test_throughput_concurrency_is_cores_over_threads(self):
        runs = plan_runs([8, 12, 24], cores=24, modes=["throughput"], reps=1)
        conc = {r["threads"]: r["concurrency"] for r in runs}
        assert conc == {8: 3, 12: 2, 24: 1}

    def test_single_mode_is_always_one_job(self):
        runs = plan_runs([8, 12, 24], cores=24, modes=["single"], reps=1)
        assert all(r["concurrency"] == 1 for r in runs)

    def test_reps_and_both_modes_multiply(self):
        runs = plan_runs([8, 16], cores=32, modes=["single", "throughput"], reps=3)
        assert len(runs) == 2 * 2 * 3  # modes x threads x reps
        assert {r["rep"] for r in runs} == {1, 2, 3}

    def test_concurrency_never_zero(self):
        # threads > cores -> still one job, not zero
        runs = plan_runs([64], cores=24, modes=["throughput"], reps=1)
        assert runs[0]["concurrency"] == 1


class TestSummarize:
    def _rows(self, data):
        # data: list of (mode, threads, throughput)
        return [
            {"mode": m, "threads": t, "throughput_reads_s": tp}
            for (m, t, tp) in data
        ]

    def test_picks_highest_throughput_split(self):
        rows = self._rows([
            ("throughput", 24, 100), ("throughput", 24, 110),  # median 105
            ("throughput", 12, 200), ("throughput", 12, 190),  # median 195  <- best
            ("throughput", 8, 150), ("throughput", 8, 150),    # median 150
        ])
        s = summarize(rows)
        assert s["best_threads"] == 12
        assert s["throughput_by_threads"][12] == 195

    def test_single_mode_excluded_from_best(self):
        rows = self._rows([
            ("single", 24, 9999),          # single must not win
            ("throughput", 12, 100),
            ("throughput", 8, 120),
        ])
        s = summarize(rows)
        assert s["best_threads"] == 8

    def test_empty_is_safe(self):
        assert summarize([])["best_threads"] is None
