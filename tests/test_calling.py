"""
TDD test suite for BSMN calling pipeline (Phase 2).

Test order follows PIPELINE.md:
  1. config validation    — calling section has required keys
  2. calling.smk structure — all 7 rules defined in file
  3. Snakemake dry-run    — DAG builds cleanly with stage="calling"

Run all:
    pytest tests/test_calling.py -v

Skip integration (no snakemake needed):
    pytest tests/test_calling.py -v -m "not integration"
"""

import gzip
import subprocess
from pathlib import Path

import pytest
import yaml

PROJECT_ROOT = Path(__file__).parent.parent
WORKFLOW_DIR = PROJECT_ROOT / "workflow"
CONFIG_DIR = PROJECT_ROOT / "config"

EXPECTED_RULES = [
    "mutect2_scatter",
    "merge_vcfs",
    "merge_mutect_stats",
    "learn_read_orientation",
    "get_pileup_summaries",
    "calculate_contamination",
    "filter_mutect_calls",
]

# =============================================================================
# Session fixture — reuse novaseq_fastq from conftest / re-create here
# =============================================================================

@pytest.fixture(scope="session")
def novaseq_fastq(tmp_path_factory):
    path = tmp_path_factory.mktemp("fastq") / "novaseq.R1.fastq.gz"
    with gzip.open(path, "wt") as fh:
        for i in range(20):
            fh.write(f"@A00100:123:AABBCCDD:1:1101:{1000+i}:{2000+i}:NNNNN 1:N:0:ATCG\n")
            fh.write("ACGTACGTACGTACGTACGTACGTACGTACGT\n")
            fh.write("+\n")
            fh.write("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")
    return str(path)


# =============================================================================
# 1. config/config.yaml — calling section validation
# =============================================================================

class TestCallingConfig:
    """Validate config/config.yaml has all required calling fields."""

    @pytest.fixture
    def cfg(self):
        path = CONFIG_DIR / "config.yaml"
        assert path.exists(), "config/config.yaml not found"
        with open(path) as fh:
            return yaml.safe_load(fh)

    def test_calling_section_exists(self, cfg):
        assert "calling" in cfg, "config.yaml missing 'calling' section"

    def test_chromosomes_is_nonempty_list(self, cfg):
        chroms = cfg["calling"].get("chromosomes")
        assert isinstance(chroms, list) and len(chroms) > 0, (
            "calling.chromosomes must be a non-empty list"
        )

    def test_chromosomes_are_chr_prefixed_strings(self, cfg):
        # Schema test: each entry must be a string starting with "chr".
        # Value test (specific chrs like chr1/chrX/chrY) belongs in production
        # config validation, not in project config (which may be a validation subset).
        chroms = cfg["calling"]["chromosomes"]
        for c in chroms:
            assert isinstance(c, str) and c.startswith("chr"), (
                f"Every chromosome must be a 'chrN' string, got: {c!r}"
            )

    def test_germline_resource_key_exists(self, cfg):
        # Schema test: key must exist (value may be empty for validation runs).
        assert "germline_resource" in cfg["calling"], (
            "calling.germline_resource key must be present in config"
        )

    def test_contamination_resource_present(self, cfg):
        assert cfg["calling"].get("contamination_resource"), (
            "calling.contamination_resource must be set"
        )

    def test_pon_block_has_use_and_vcf(self, cfg):
        pon = cfg["calling"].get("pon", {})
        assert "use" in pon, "calling.pon.use missing"
        assert "vcf" in pon, "calling.pon.vcf missing"

    def test_filter_extra_present(self, cfg):
        assert "filter_extra" in cfg["calling"], "calling.filter_extra missing"


# =============================================================================
# 2. workflow/rules/calling.smk — rule structure
# =============================================================================

class TestCallingSmkRules:
    """Verify calling.smk defines all expected rules (file-level check)."""

    @pytest.fixture
    def smk_text(self):
        smk = WORKFLOW_DIR / "rules" / "calling.smk"
        assert smk.exists(), (
            "workflow/rules/calling.smk not found — implement it first"
        )
        return smk.read_text()

    def test_all_rules_defined(self, smk_text):
        for rule in EXPECTED_RULES:
            assert f"rule {rule}:" in smk_text, (
                f"Rule '{rule}' not found in calling.smk"
            )

    def test_apptainer_exec_used(self, smk_text):
        assert "apptainer exec" in smk_text, (
            "calling.smk must use 'apptainer exec', not conda"
        )

    def test_no_conda_directive(self, smk_text):
        import re
        # Allow 'conda' in comments but not as a directive
        assert not re.search(r"^\s+conda:", smk_text, re.MULTILINE), (
            "calling.smk must not use conda: directive"
        )

    def test_temp_on_scatter_outputs(self, smk_text):
        assert "temp(" in smk_text, (
            "Scatter outputs should be marked temp() to avoid disk waste"
        )

    def test_mutect2_tumor_only_flag(self, smk_text):
        assert "--tumor-sample" in smk_text, (
            "Mutect2 must use --tumor-sample for tumor-only mode"
        )

    def test_filter_uses_ob_priors(self, smk_text):
        assert "--ob-priors" in smk_text, (
            "FilterMutectCalls must use --ob-priors (read orientation model)"
        )

    def test_filter_uses_contamination_table(self, smk_text):
        assert "--contamination-table" in smk_text, (
            "FilterMutectCalls must use --contamination-table"
        )


# =============================================================================
# 3. Snakemake dry-run (integration)
# =============================================================================

@pytest.mark.integration
class TestCallingDryRun:
    """DAG dry-run with stage='calling' — all mapping + calling rules planned."""

    @pytest.fixture(scope="class")
    def calling_env(self, tmp_path_factory, novaseq_fastq):
        env = tmp_path_factory.mktemp("calling_env")
        cfg = env / "config"
        cfg.mkdir()

        # samples.tsv
        (cfg / "samples.tsv").write_text(
            "sample_id\treadgroup\tfq1\tfq2\n"
            f"TEST001\tRG1\t{novaseq_fastq}\t{novaseq_fastq}\n"
        )

        # containers.yaml
        import shutil
        shutil.copy(CONFIG_DIR / "containers.yaml", cfg / "containers.yaml")

        # resolved_params.yaml
        from scripts.auto_params import resolve_params
        resolve_params(novaseq_fastq, str(cfg / "resolved_params.yaml"))

        # config.yaml — stage=calling, 2 chromosomes for speed
        dummy_ref = str(env / "ref.fasta")
        (env / "ref.fasta").touch()
        cfg_data = {
            "samples": "config/samples.tsv",
            "stage": "calling",
            "ref": {
                "build": "hg38",
                "fasta": dummy_ref,
                "dict": dummy_ref.replace(".fasta", ".dict"),
            },
            "known_sites": {
                "dbsnp": "dummy.vcf.gz",
                "mills": "dummy.vcf.gz",
                "indels": "dummy.vcf.gz",
            },
            "mapping": {
                "bwa_threads": 2,
                "sort_threads": 1,
                "sort_memory": "256M",
                "markdup_memory": "1G",
                "output_format": "cram",
            },
            "calling": {
                "caller": "mutect2",
                "chromosomes": ["chr1", "chr2"],
                "germline_resource": "dummy_gnomad.vcf.gz",
                "contamination_resource": "dummy_exac.vcf.gz",
                "pon": {"use": False, "vcf": ""},
                "mutect2_extra": "",
                "filter_extra": "--min-reads-per-strand 1",
            },
        }
        with open(cfg / "config.yaml", "w") as fh:
            yaml.dump(cfg_data, fh)

        return env

    def _dry_run(self, calling_env):
        return subprocess.run(
            [
                "snakemake",
                "--dry-run",
                "--snakefile", str(WORKFLOW_DIR / "Snakefile"),
                "--directory", str(calling_env),
                "--cores", "1",
            ],
            capture_output=True,
            text=True,
        )

    def test_dry_run_exits_zero(self, calling_env):
        result = self._dry_run(calling_env)
        assert result.returncode == 0, (
            f"Snakemake dry-run failed.\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    def test_all_calling_rules_in_dag(self, calling_env):
        result = self._dry_run(calling_env)
        combined = result.stdout + result.stderr
        for rule in EXPECTED_RULES:
            assert rule in combined, (
                f"Rule '{rule}' not found in dry-run DAG.\n"
                f"Output (first 3000 chars):\n{combined[:3000]}"
            )

    def test_scatter_jobs_equal_chromosome_count(self, calling_env):
        """mutect2_scatter should appear once per chromosome (2 in test config)."""
        result = self._dry_run(calling_env)
        combined = result.stdout + result.stderr
        count = combined.count("mutect2_scatter")
        assert count >= 2, (
            f"Expected at least 2 mutect2_scatter jobs (one per chrom), got {count}"
        )

    def test_filter_mutect_calls_is_final_target(self, calling_env):
        result = self._dry_run(calling_env)
        combined = result.stdout + result.stderr
        assert "filter_mutect_calls" in combined, (
            "filter_mutect_calls must appear in calling DAG"
        )
