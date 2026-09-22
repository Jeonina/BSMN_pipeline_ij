"""Portability of run.py's host-specific wiring.

Two settings decide whether a clone runs correctly on someone else's server, and
both fail quietly when absent: the Snakemake memory budget (without it every
`resources: mem_mb` declaration is inert) and the Apptainer bind list (without it
a dry-run passes and the real run dies inside the container).
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run  # noqa: E402


def _covered(binds: list[str], target: Path) -> bool:
    """True if some bind exposes `target` (itself or an ancestor)."""
    target = target.resolve()
    return any(target == Path(b) or Path(b) in target.parents for b in binds)


@pytest.fixture
def workdir(tmp_path, monkeypatch):
    """A working directory with no ambient binds, so assertions are exact."""
    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.chdir(work)
    monkeypatch.setattr(run, "_ALWAYS_BIND", ())
    return work


class TestApptainerBinds:
    def test_input_on_another_mount_is_bound(self, workdir, tmp_path):
        data = tmp_path / "storage" / "wgs"
        data.mkdir(parents=True)
        rows = [{"fq1": str(data / "s_R1.fq.gz"), "fq2": str(data / "s_R2.fq.gz")}]

        assert _covered(run._apptainer_binds(rows), data)

    def test_input_inside_workdir_adds_nothing(self, workdir):
        inside = workdir / "fastq"
        inside.mkdir()
        rows = [{"fq1": str(inside / "s_R1.fq.gz")}]

        # Already visible to the container; must not enlarge the bind list.
        assert run._apptainer_binds(rows) == run._apptainer_binds([])

    def test_results_symlink_target_is_bound(self, workdir, tmp_path):
        nfs = tmp_path / "nfs" / "results"
        nfs.mkdir(parents=True)
        (workdir / "results").symlink_to(nfs)

        # The production layout: outputs land on NFS through a symlink, so the
        # symlink TARGET is what the container needs, not the link.
        assert _covered(run._apptainer_binds([]), nfs)

    def test_bam_mode_rows_are_bound(self, workdir, tmp_path):
        data = tmp_path / "crams"
        data.mkdir()
        rows = [{"bam": str(data / "S1.cram")}]

        assert _covered(run._apptainer_binds(rows), data)

    def test_nested_paths_collapse_to_the_ancestor(self, workdir, tmp_path):
        parent = tmp_path / "storage"
        child = parent / "run1"
        child.mkdir(parents=True)
        rows = [{"fq1": str(parent / "a.fq.gz"), "fq2": str(child / "b.fq.gz")}]

        binds = run._apptainer_binds(rows)
        assert str(child) not in binds
        assert _covered(binds, child)

    def test_missing_input_paths_do_not_raise(self, workdir, tmp_path):
        rows = [{"fq1": str(tmp_path / "gone" / "x.fq.gz")}, {}, {"fq1": None}]

        run._apptainer_binds(rows)  # must not raise


class TestMemoryBudget:
    def test_explicit_value_wins(self):
        assert run._resolve_mem_mb(40000) == 40000

    def test_derived_from_host_when_unset(self):
        budget = run._resolve_mem_mb(None)

        assert isinstance(budget, int)
        assert budget >= 4096
