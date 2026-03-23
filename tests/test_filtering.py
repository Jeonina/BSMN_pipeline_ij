"""
TDD test suite for BSMN filtering pipeline (Phase 3).

Test order follows PIPELINE.md:
  1. config validation      — filtering section has required keys
  2. filtering.smk structure — all 4 rules defined in file
  3. Script unit tests       — accessibility, germline, VAF, PON logic
  4. Snakemake dry-run       — DAG builds cleanly with stage="filtering"

Run all:
    pytest tests/test_filtering.py -v

Skip integration (no snakemake needed):
    pytest tests/test_filtering.py -v -m "not integration"
"""

import gzip
import io
import subprocess
from pathlib import Path

import pytest
import yaml

PROJECT_ROOT = Path(__file__).parent.parent
WORKFLOW_DIR = PROJECT_ROOT / "workflow"
CONFIG_DIR = PROJECT_ROOT / "config"
SCRIPTS_DIR = PROJECT_ROOT / "scripts"

EXPECTED_FILTER_RULES = [
    "accessibility_filter",
    "germline_filter",
    "vaf_filter",
    "pon_mask_filter",
]

# =============================================================================
# Session fixture — reuse novaseq_fastq from conftest
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
# 1. config/config.yaml — filtering section validation
# =============================================================================


class TestFilteringConfig:
    """Validate config/config.yaml has all required filtering fields."""

    @pytest.fixture
    def cfg(self):
        path = CONFIG_DIR / "config.yaml"
        assert path.exists(), "config/config.yaml not found"
        with open(path) as fh:
            return yaml.safe_load(fh)

    def test_filtering_section_exists(self, cfg):
        assert "filtering" in cfg, "config.yaml missing 'filtering' section"

    def test_gnomad_snps_key_present(self, cfg):
        gnomad = cfg["filtering"].get("gnomad", {})
        assert "snps" in gnomad, "filtering.gnomad.snps missing"

    def test_gnomad_af_threshold_present(self, cfg):
        gnomad = cfg["filtering"].get("gnomad", {})
        assert "af_threshold" in gnomad, "filtering.gnomad.af_threshold missing"

    def test_vaf_block_present(self, cfg):
        vaf = cfg["filtering"].get("vaf", {})
        assert "min_alt_count" in vaf, "filtering.vaf.min_alt_count missing"
        assert "p_binom_threshold" in vaf, "filtering.vaf.p_binom_threshold missing"

    def test_pon_mask_fasta_present(self, cfg):
        pon = cfg["filtering"].get("pon_mask", {})
        assert "fasta" in pon, "filtering.pon_mask.fasta missing"

    def test_mask_1kg_present(self, cfg):
        assert cfg["filtering"].get("mask_1kg"), "filtering.mask_1kg must be set"


# =============================================================================
# 2. workflow/rules/filtering.smk — rule structure
# =============================================================================


class TestFilteringSmkRules:
    """Verify filtering.smk defines all expected rules."""

    @pytest.fixture
    def smk_text(self):
        smk = WORKFLOW_DIR / "rules" / "filtering.smk"
        assert smk.exists(), (
            "workflow/rules/filtering.smk not found — implement it first"
        )
        return smk.read_text()

    def test_all_rules_defined(self, smk_text):
        for rule in EXPECTED_FILTER_RULES:
            assert f"rule {rule}:" in smk_text, (
                f"Rule '{rule}' not found in filtering.smk"
            )

    def test_apptainer_exec_used(self, smk_text):
        assert "apptainer exec" in smk_text, (
            "filtering.smk must use 'apptainer exec', not conda"
        )

    def test_no_conda_directive(self, smk_text):
        import re
        assert not re.search(r"^\s+conda:", smk_text, re.MULTILINE), (
            "filtering.smk must not use conda: directive"
        )

    def test_temp_on_intermediate_outputs(self, smk_text):
        assert "temp(" in smk_text, (
            "Intermediate filter outputs should be marked temp() to avoid disk waste"
        )

    def test_samtools_faidx_in_accessibility_rule(self, smk_text):
        assert "samtools faidx" in smk_text, (
            "accessibility_filter must call samtools faidx for mask lookup"
        )

    def test_vaf_filter_has_min_mapq(self, smk_text):
        assert "--min-mapq" in smk_text, (
            "vaf_filter must pass --min-mapq (mapping quality threshold)"
        )

    def test_pon_mask_uses_pon_fasta(self, smk_text):
        assert "pon_fasta" in smk_text or "pon-fasta" in smk_text, (
            "pon_mask_filter must reference pon_fasta resource"
        )

    def test_final_output_not_temp(self, smk_text):
        import re
        # final.txt should NOT be wrapped in temp()
        final_line = [l for l in smk_text.splitlines() if "final.txt" in l]
        assert final_line, "filtering.smk must define final.txt output"
        for line in final_line:
            assert "temp(" not in line, "final.txt output must NOT be temp()"


# =============================================================================
# 3a. scripts/germline_filter.py — unit tests
# =============================================================================


class TestGermlineFilterScript:
    """Unit tests for scripts/germline_filter.py."""

    @pytest.fixture
    def module(self):
        import sys
        sys.path.insert(0, str(PROJECT_ROOT))
        from scripts.germline_filter import (
            load_germline_set,
            is_germline,
            filter_txt,
        )
        return load_germline_set, is_germline, filter_txt

    def _make_gnomad_gz(self, tmp_path, entries):
        """Write a gzipped gnomAD SNP file (chrom\tpos\tref\talt, no chr prefix)."""
        gz_path = tmp_path / "gnomad.txt.gz"
        with gzip.open(gz_path, "wt") as fh:
            for entry in entries:
                fh.write("\t".join(entry) + "\n")
        return str(gz_path)

    def test_load_germline_set_basic(self, module, tmp_path):
        load_germline_set, _, _ = module
        gz = self._make_gnomad_gz(tmp_path, [
            ["1", "12345", "A", "G"],
            ["X", "67890", "C", "T"],
        ])
        known = load_germline_set(gz)
        assert "1:12345:A:G" in known
        assert "X:67890:C:T" in known

    def test_is_germline_returns_true_for_known(self, module, tmp_path):
        load_germline_set, is_germline, _ = module
        gz = self._make_gnomad_gz(tmp_path, [["1", "12345", "A", "G"]])
        known = load_germline_set(gz)
        assert is_germline("chr1", "12345", "A", "G", known) is True

    def test_is_germline_strips_chr_prefix(self, module, tmp_path):
        load_germline_set, is_germline, _ = module
        gz = self._make_gnomad_gz(tmp_path, [["1", "100", "T", "C"]])
        known = load_germline_set(gz)
        # Both with and without chr prefix should match
        assert is_germline("chr1", "100", "T", "C", known) is True
        assert is_germline("1", "100", "T", "C", known) is True

    def test_is_germline_returns_false_for_unknown(self, module, tmp_path):
        load_germline_set, is_germline, _ = module
        gz = self._make_gnomad_gz(tmp_path, [["1", "12345", "A", "G"]])
        known = load_germline_set(gz)
        assert is_germline("chr1", "99999", "A", "G", known) is False

    def test_filter_txt_removes_germline_variants(self, module, tmp_path):
        load_germline_set, _, filter_txt = module
        gz = self._make_gnomad_gz(tmp_path, [["1", "12345", "A", "G"]])
        known = load_germline_set(gz)

        lines = (
            "chr1\t12345\tA\tG\n"   # germline → removed
            "chr1\t99999\tC\tT\n"   # somatic → kept
        )
        infile = io.StringIO(lines)
        outfile = io.StringIO()
        filter_txt(infile, known, outfile)
        result = outfile.getvalue()
        assert "12345" not in result
        assert "99999" in result

    def test_filter_txt_preserves_comment_lines(self, module, tmp_path):
        load_germline_set, _, filter_txt = module
        gz = self._make_gnomad_gz(tmp_path, [])
        known = load_germline_set(gz)

        lines = "#header line\nchr1\t100\tA\tG\n"
        infile = io.StringIO(lines)
        outfile = io.StringIO()
        filter_txt(infile, known, outfile)
        assert "#header line" in outfile.getvalue()


# =============================================================================
# 3b. scripts/vaf_filter.py — unit tests (pure function only)
# =============================================================================


class TestVafFilterScript:
    """Unit tests for scripts/vaf_filter.py — pure function passes_vaf_filter."""

    @pytest.fixture
    def passes_vaf_filter(self):
        import sys
        sys.path.insert(0, str(PROJECT_ROOT))
        from scripts.vaf_filter import passes_vaf_filter
        return passes_vaf_filter

    def test_somatic_mosaic_passes(self, passes_vaf_filter):
        """Low VAF (10/100 = 0.1) with >=5 alt reads should pass."""
        assert passes_vaf_filter(alt_n=10, depth=100) is True

    def test_germline_hetero_fails(self, passes_vaf_filter):
        """VAF ~0.5 (50/100) is likely germline → should fail."""
        assert passes_vaf_filter(alt_n=50, depth=100) is False

    def test_low_alt_count_fails(self, passes_vaf_filter):
        """Even with good p-value, alt_n < min_alt → fails."""
        assert passes_vaf_filter(alt_n=3, depth=1000) is False

    def test_zero_alt_fails(self, passes_vaf_filter):
        assert passes_vaf_filter(alt_n=0, depth=100) is False

    def test_zero_depth_fails(self, passes_vaf_filter):
        assert passes_vaf_filter(alt_n=0, depth=0) is False

    def test_very_low_vaf_passes(self, passes_vaf_filter):
        """Very low mosaic (5/300 ≈ 1.7%) with exactly min_alt reads."""
        assert passes_vaf_filter(alt_n=5, depth=300) is True

    def test_custom_thresholds(self, passes_vaf_filter):
        """Custom p_threshold and min_alt are respected."""
        # min_alt=10: alt_n=7 should fail
        assert passes_vaf_filter(alt_n=7, depth=200, min_alt=10) is False
        # p_threshold=0.5: should be much more lenient
        assert passes_vaf_filter(alt_n=10, depth=30, p_threshold=0.5) is True


# =============================================================================
# 3c. scripts/pon_mask_filter.py — unit tests
# =============================================================================


class TestPonMaskFilterScript:
    """Unit tests for scripts/pon_mask_filter.py — IUPAC logic."""

    @pytest.fixture
    def pon_funcs(self):
        import sys
        sys.path.insert(0, str(PROJECT_ROOT))
        from scripts.pon_mask_filter import iupac_match, pon_passes
        return iupac_match, pon_passes

    # --- Exact base matches (single nucleotide) ---
    def test_exact_match_A_fails(self, pon_funcs):
        _, pon_passes = pon_funcs
        assert pon_passes("A", "A") is False

    def test_exact_mismatch_passes(self, pon_funcs):
        _, pon_passes = pon_funcs
        assert pon_passes("G", "A") is True

    # --- Degenerate IUPAC codes ---
    def test_R_matches_A_and_G(self, pon_funcs):
        iupac_match, _ = pon_funcs
        assert iupac_match("A", "R") is True   # R = A|G
        assert iupac_match("G", "R") is True
        assert iupac_match("C", "R") is False

    def test_Y_matches_C_and_T(self, pon_funcs):
        iupac_match, _ = pon_funcs
        assert iupac_match("C", "Y") is True   # Y = C|T
        assert iupac_match("T", "Y") is True
        assert iupac_match("A", "Y") is False

    def test_S_matches_G_and_C(self, pon_funcs):
        iupac_match, _ = pon_funcs
        assert iupac_match("G", "S") is True
        assert iupac_match("C", "S") is True
        assert iupac_match("A", "S") is False

    def test_W_matches_A_and_T(self, pon_funcs):
        iupac_match, _ = pon_funcs
        assert iupac_match("A", "W") is True
        assert iupac_match("T", "W") is True
        assert iupac_match("G", "W") is False

    def test_K_matches_G_and_T(self, pon_funcs):
        iupac_match, _ = pon_funcs
        assert iupac_match("G", "K") is True
        assert iupac_match("T", "K") is True
        assert iupac_match("A", "K") is False

    def test_M_matches_A_and_C(self, pon_funcs):
        iupac_match, _ = pon_funcs
        assert iupac_match("A", "M") is True
        assert iupac_match("C", "M") is True
        assert iupac_match("G", "M") is False

    def test_B_matches_C_G_T(self, pon_funcs):
        iupac_match, _ = pon_funcs
        assert iupac_match("C", "B") is True
        assert iupac_match("G", "B") is True
        assert iupac_match("T", "B") is True
        assert iupac_match("A", "B") is False

    def test_D_matches_A_G_T(self, pon_funcs):
        iupac_match, _ = pon_funcs
        assert iupac_match("A", "D") is True
        assert iupac_match("G", "D") is True
        assert iupac_match("T", "D") is True
        assert iupac_match("C", "D") is False

    def test_H_matches_A_C_T(self, pon_funcs):
        iupac_match, _ = pon_funcs
        assert iupac_match("A", "H") is True
        assert iupac_match("C", "H") is True
        assert iupac_match("T", "H") is True
        assert iupac_match("G", "H") is False

    def test_V_matches_A_C_G(self, pon_funcs):
        iupac_match, _ = pon_funcs
        assert iupac_match("A", "V") is True
        assert iupac_match("C", "V") is True
        assert iupac_match("G", "V") is True
        assert iupac_match("T", "V") is False

    # --- Special cases ---
    def test_N_always_fails(self, pon_funcs):
        _, pon_passes = pon_funcs
        for base in ("A", "C", "G", "T"):
            assert pon_passes(base, "N") is False, f"N should always fail for alt={base}"

    def test_star_always_passes(self, pon_funcs):
        _, pon_passes = pon_funcs
        for base in ("A", "C", "G", "T"):
            assert pon_passes(base, "*") is True, f"* should always pass for alt={base}"

    def test_unknown_iupac_base_is_conservative_fail(self, pon_funcs):
        iupac_match, _ = pon_funcs
        # Unknown/unexpected base: conservative default is Fail (returns True from iupac_match)
        assert iupac_match("A", "?") is True

    def test_case_insensitive(self, pon_funcs):
        _, pon_passes = pon_funcs
        # Lowercase alt
        assert pon_passes("a", "A") is False
        assert pon_passes("g", "R") is False


# =============================================================================
# 4. Snakemake dry-run with stage="filtering" (integration)
# =============================================================================


@pytest.mark.integration
class TestFilteringDryRun:
    """DAG dry-run with stage='filtering' — all rules planned."""

    @pytest.fixture(scope="class")
    def filtering_env(self, tmp_path_factory, novaseq_fastq):
        env = tmp_path_factory.mktemp("filtering_env")
        cfg_dir = env / "config"
        cfg_dir.mkdir()

        # samples.tsv
        (cfg_dir / "samples.tsv").write_text(
            "sample_id\treadgroup\tfq1\tfq2\n"
            f"TEST001\tRG1\t{novaseq_fastq}\t{novaseq_fastq}\n"
        )

        # containers.yaml
        import shutil
        shutil.copy(CONFIG_DIR / "containers.yaml", cfg_dir / "containers.yaml")

        # resolved_params.yaml
        from scripts.auto_params import resolve_params
        resolve_params(novaseq_fastq, str(cfg_dir / "resolved_params.yaml"))

        # config.yaml — stage=filtering, minimal chromosomes
        dummy_ref = str(env / "ref.fasta")
        (env / "ref.fasta").touch()
        cfg_data = {
            "samples": "config/samples.tsv",
            "stage": "filtering",
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
            "filtering": {
                "gnomad": {
                    "snps": "dummy_gnomad_snps.txt.gz",
                    "af_threshold": 0.001,
                },
                "vaf": {
                    "min_alt_count": 5,
                    "p_binom_threshold": 1.0e-6,
                },
                "pon_mask": {"fasta": "dummy_pon.fa"},
                "mask_1kg": "dummy_1kg_mask.fa.gz",
            },
        }
        with open(cfg_dir / "config.yaml", "w") as fh:
            yaml.dump(cfg_data, fh)

        return env

    def _dry_run(self, filtering_env):
        return subprocess.run(
            [
                "snakemake",
                "--dry-run",
                "--snakefile", str(WORKFLOW_DIR / "Snakefile"),
                "--directory", str(filtering_env),
                "--cores", "1",
            ],
            capture_output=True,
            text=True,
        )

    def test_dry_run_exits_zero(self, filtering_env):
        result = self._dry_run(filtering_env)
        assert result.returncode == 0, (
            f"Snakemake dry-run failed.\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    def test_all_filtering_rules_in_dag(self, filtering_env):
        result = self._dry_run(filtering_env)
        combined = result.stdout + result.stderr
        for rule in EXPECTED_FILTER_RULES:
            assert rule in combined, (
                f"Rule '{rule}' not found in dry-run DAG.\n"
                f"Output (first 3000 chars):\n{combined[:3000]}"
            )

    def test_pon_mask_filter_is_final_target(self, filtering_env):
        result = self._dry_run(filtering_env)
        combined = result.stdout + result.stderr
        assert "pon_mask_filter" in combined, (
            "pon_mask_filter must appear as final filtering step in DAG"
        )

    def test_calling_rules_also_in_dag(self, filtering_env):
        """Filtering depends on calling, which depends on mapping — all in DAG."""
        result = self._dry_run(filtering_env)
        combined = result.stdout + result.stderr
        assert "mutect2_scatter" in combined, (
            "mutect2_scatter (from calling stage) must appear in filtering DAG"
        )
