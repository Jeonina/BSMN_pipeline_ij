"""
TDD test suite for BSMN mapping pipeline (Phase 1).

Test order follows PIPELINE.md spec:
  1. auto_params unit tests   — RED first, then GREEN after scripts/auto_params.py exists
  2. containers.yaml validation
  3. Snakemake dry-run (DAG check)  — marked @pytest.mark.integration

Run all:
    pytest tests/test_mapping.py -v

Skip integration tests (no snakemake install required):
    pytest tests/test_mapping.py -v -m "not integration"
"""

import gzip
import os
import subprocess
import tempfile
from pathlib import Path

import pytest
import yaml

PROJECT_ROOT = Path(__file__).parent.parent
WORKFLOW_DIR = PROJECT_ROOT / "workflow"
CONFIG_DIR = PROJECT_ROOT / "config"
TESTS_DIR = PROJECT_ROOT / "tests"
DATA_DIR = TESTS_DIR / "data"


# =============================================================================
# Session fixtures — synthetic FASTQ data
# =============================================================================

@pytest.fixture(scope="session")
def novaseq_fastq(tmp_path_factory):
    """
    Tiny synthetic NovaSeq R1 FASTQ.
    Read name format: @instrument:run:flowcell:lane:tile:x:y (7 colons → NovaSeq).
    """
    path = tmp_path_factory.mktemp("fastq") / "novaseq.R1.fastq.gz"
    with gzip.open(path, "wt") as fh:
        for i in range(20):
            # NovaSeq with UMI: 7 colons in read name (instrument:run:flowcell:lane:tile:x:y:UMI)
            fh.write(f"@A00100:123:AABBCCDD:1:1101:{1000 + i}:{2000 + i}:NNNNN 1:N:0:ATCG\n")
            fh.write("ACGTACGTACGTACGTACGTACGTACGTACGT\n")
            fh.write("+\n")
            fh.write("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")
    return str(path)


@pytest.fixture(scope="session")
def hiseq_fastq(tmp_path_factory):
    """
    Tiny synthetic HiSeq R1 FASTQ.
    Read name format: @instrument:lane:tile:x:y (4 colons → HiSeq).
    """
    path = tmp_path_factory.mktemp("fastq") / "hiseq.R1.fastq.gz"
    with gzip.open(path, "wt") as fh:
        for i in range(20):
            fh.write(f"@HISEQ2500:1:1101:{1000 + i}:{2000 + i} 1:N:0:\n")
            fh.write("ACGTACGTACGTACGTACGTACGTACGTACGT\n")
            fh.write("+\n")
            fh.write("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")
    return str(path)


# =============================================================================
# 1. auto_params — sequencer detection
# =============================================================================

class TestDetectSequencer:
    """RED: written before auto_params.py exists."""

    def test_novaseq_7_colons_returns_novaseq(self, novaseq_fastq):
        from scripts.auto_params import detect_sequencer
        assert detect_sequencer(novaseq_fastq) == "NovaSeq"

    def test_hiseq_fewer_colons_returns_hiseq(self, hiseq_fastq):
        from scripts.auto_params import detect_sequencer
        assert detect_sequencer(hiseq_fastq) == "HiSeq"

    def test_uncompressed_fastq_works(self, tmp_path):
        """detect_sequencer must handle plain (non-gz) FASTQ."""
        from scripts.auto_params import detect_sequencer
        fq = tmp_path / "plain.fastq"
        fq.write_text(
            "@A00100:123:AABBCCDD:1:1101:1000:2000:NNNNN 1:N:0:ATCG\nACGT\n+\nIIII\n"
        )
        assert detect_sequencer(str(fq)) == "NovaSeq"

    def test_boundary_exactly_7_colons_is_novaseq(self, tmp_path):
        """Exactly 7 colons is the NovaSeq boundary."""
        from scripts.auto_params import detect_sequencer
        fq = tmp_path / "boundary.fastq.gz"
        # 7 colons: a:b:c:d:e:f:g:h → 7 colons
        with gzip.open(fq, "wt") as f:
            f.write("@a:b:c:d:e:f:g:h 1:N:0:\nACGT\n+\nIIII\n")
        assert detect_sequencer(str(fq)) == "NovaSeq"

    def test_boundary_6_colons_is_hiseq(self, tmp_path):
        """6 colons → HiSeq."""
        from scripts.auto_params import detect_sequencer
        fq = tmp_path / "hiseq6.fastq.gz"
        with gzip.open(fq, "wt") as f:
            f.write("@a:b:c:d:e:f:g 1:N:0:\nACGT\n+\nIIII\n")
        assert detect_sequencer(str(fq)) == "HiSeq"


# =============================================================================
# 2. auto_params — ODPD selection
# =============================================================================

class TestGetODPD:
    def test_novaseq_odpd_is_2500(self):
        from scripts.auto_params import get_optical_duplicate_pixel_distance
        assert get_optical_duplicate_pixel_distance("NovaSeq") == 2500

    def test_hiseq_odpd_is_100(self):
        from scripts.auto_params import get_optical_duplicate_pixel_distance
        assert get_optical_duplicate_pixel_distance("HiSeq") == 100

    def test_unknown_sequencer_defaults_to_hiseq(self):
        from scripts.auto_params import get_optical_duplicate_pixel_distance
        assert get_optical_duplicate_pixel_distance("MiSeq") == 100


# =============================================================================
# 3. auto_params — system resource queries
# =============================================================================

class TestGetBwaThreads:
    def test_returns_positive_int(self):
        from scripts.auto_params import get_bwa_threads
        t = get_bwa_threads()
        assert isinstance(t, int)
        assert t >= 1


class TestGetBqsrMemory:
    def test_returns_positive_int(self):
        from scripts.auto_params import get_bqsr_memory_gb
        m = get_bqsr_memory_gb()
        assert isinstance(m, int)
        assert m >= 4

    def test_bounded_between_4_and_64(self):
        from scripts.auto_params import get_bqsr_memory_gb
        m = get_bqsr_memory_gb()
        assert 4 <= m <= 64


# =============================================================================
# 4. auto_params — resolve_params() end-to-end
# =============================================================================

class TestResolveParams:
    def test_writes_yaml_file(self, novaseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        out = str(tmp_path / "resolved.yaml")
        resolve_params(novaseq_fastq, out)
        assert os.path.exists(out)

    def test_yaml_is_parseable(self, novaseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        out = str(tmp_path / "resolved.yaml")
        resolve_params(novaseq_fastq, out)
        with open(out) as fh:
            data = yaml.safe_load(fh)
        assert isinstance(data, dict)

    def test_required_keys_present(self, novaseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        params = resolve_params(novaseq_fastq, str(tmp_path / "r.yaml"))
        required = [
            "resolved_at",
            "input_fastq",
            "sequencer",
            "bwa_threads",
            "sort_threads",
            "sort_memory",
            "optical_duplicate_pixel_distance",
            "bqsr_memory_gb",
            "markdup_memory",
        ]
        for key in required:
            assert key in params, f"Missing key in resolved_params: {key}"

    def test_novaseq_sets_odpd_2500(self, novaseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        params = resolve_params(novaseq_fastq, str(tmp_path / "r.yaml"))
        assert params["optical_duplicate_pixel_distance"] == 2500

    def test_hiseq_sets_odpd_100(self, hiseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        params = resolve_params(hiseq_fastq, str(tmp_path / "r.yaml"))
        assert params["optical_duplicate_pixel_distance"] == 100

    def test_creates_parent_directory(self, novaseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        out = str(tmp_path / "nested" / "subdir" / "resolved.yaml")
        resolve_params(novaseq_fastq, out)
        assert os.path.exists(out)

    def test_input_fastq_stored_as_absolute_path(self, novaseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        params = resolve_params(novaseq_fastq, str(tmp_path / "r.yaml"))
        assert os.path.isabs(params["input_fastq"])

    def test_resolved_at_is_iso_format(self, novaseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        import datetime
        params = resolve_params(novaseq_fastq, str(tmp_path / "r.yaml"))
        # Should parse without raising
        datetime.datetime.fromisoformat(params["resolved_at"])


# =============================================================================
# 5. containers.yaml structure validation
# =============================================================================

class TestContainersYaml:
    """Validate config/containers.yaml before any Snakemake run."""

    REQUIRED_TOOLS = ["bwa", "sambamba", "picard", "gatk", "samtools"]
    REQUIRED_FIELDS = ["name", "version", "uri", "sif"]

    @pytest.fixture
    def containers(self):
        path = CONFIG_DIR / "containers.yaml"
        assert path.exists(), f"config/containers.yaml not found at {path}"
        with open(path) as fh:
            return yaml.safe_load(fh)

    def test_required_tools_present(self, containers):
        for tool in self.REQUIRED_TOOLS:
            assert tool in containers, f"Missing tool in containers.yaml: {tool}"

    def test_each_tool_has_required_fields(self, containers):
        for tool in self.REQUIRED_TOOLS:
            for field in self.REQUIRED_FIELDS:
                assert field in containers[tool], (
                    f"containers.yaml[{tool}] missing field: '{field}'"
                )

    def test_sif_paths_under_containers_dir(self, containers):
        for tool, spec in containers.items():
            assert spec["sif"].startswith("containers/"), (
                f"{tool}.sif should start with 'containers/', got: {spec['sif']}"
            )

    def test_uri_starts_with_docker_or_oras(self, containers):
        for tool, spec in containers.items():
            assert spec["uri"].startswith(("docker://", "oras://", "library://")), (
                f"{tool}.uri has unexpected scheme: {spec['uri']}"
            )


# =============================================================================
# 6. Snakemake DAG dry-run  (integration — requires snakemake installed)
# =============================================================================

@pytest.mark.integration
class TestSnakemakeDryRun:
    """
    Verify the Snakemake DAG can be built without errors.
    Provides a temporary environment with minimal config + test samples.
    """

    @pytest.fixture(scope="class")
    def dry_run_env(self, tmp_path_factory, novaseq_fastq):
        """Create a temporary working directory with all required config files."""
        env = tmp_path_factory.mktemp("snakemake_env")
        cfg = env / "config"
        cfg.mkdir()

        # --- samples.tsv ---
        (cfg / "samples.tsv").write_text(
            "sample_id\treadgroup\tfq1\tfq2\n"
            f"TEST001\tRG1\t{novaseq_fastq}\t{novaseq_fastq}\n"
        )

        # --- containers.yaml (copy from project) ---
        import shutil
        shutil.copy(CONFIG_DIR / "containers.yaml", cfg / "containers.yaml")

        # --- resolved_params.yaml (generated inline) ---
        from scripts.auto_params import resolve_params
        resolve_params(novaseq_fastq, str(cfg / "resolved_params.yaml"))

        # --- minimal config.yaml ---
        dummy_ref = str(env / "ref.fasta")
        (env / "ref.fasta").touch()
        (cfg / "config.yaml").write_text(
            f"samples: config/samples.tsv\n"
            f"ref:\n"
            f"  build: hg38\n"
            f"  fasta: {dummy_ref}\n"
            f"  dict: {dummy_ref.replace('.fasta', '.dict')}\n"
            f"known_sites:\n"
            f"  dbsnp: dummy.vcf.gz\n"
            f"  mills: dummy.vcf.gz\n"
            f"  indels: dummy.vcf.gz\n"
            f"mapping:\n"
            f"  bwa_threads: 2\n"
            f"  sort_threads: 1\n"
            f"  sort_memory: 256M\n"
            f"  markdup_memory: 1G\n"
            f"  output_format: cram\n"
        )
        return env

    def test_dry_run_exits_zero(self, dry_run_env):
        result = subprocess.run(
            [
                "snakemake",
                "--dry-run",
                "--snakefile", str(WORKFLOW_DIR / "Snakefile"),
                "--directory", str(dry_run_env),
                "--cores", "1",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, (
            f"Snakemake dry-run failed.\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    def test_all_mapping_rules_present_in_dag(self, dry_run_env):
        result = subprocess.run(
            [
                "snakemake",
                "--dry-run",
                "--snakefile", str(WORKFLOW_DIR / "Snakefile"),
                "--directory", str(dry_run_env),
                "--cores", "1",
            ],
            capture_output=True,
            text=True,
        )
        combined = result.stdout + result.stderr
        expected_rules = [
            "bwa_mem_sort",
            "merge_bams",
            "mark_duplicates",
            "base_recalibrator",
            "apply_bqsr",
            "samtools_flagstat",
        ]
        for rule in expected_rules:
            assert rule in combined, (
                f"Rule '{rule}' not found in dry-run output.\n"
                f"Output:\n{combined[:2000]}"
            )
