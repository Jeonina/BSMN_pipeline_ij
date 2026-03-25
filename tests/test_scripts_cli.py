"""
TDD Phase C: CLI entry points and uncovered function tests.

Coverage targets (scripts/ directory):
  vaf_filter.py       : _clean_bases, binom_pvalue, pileup_base_counts,
                        filter_variants, main  (lines 73-199)
  pon_mask_filter.py  : query_fasta, filter_variants, main  (lines 80-141)
  germline_filter.py  : main  (lines 67-92)
  prepare_containers.py : prepare_containers, main  (lines 13-77)
  auto_params.py      : main CLI  (lines 250-290)

TDD order: tests written FIRST → confirm RED → confirm GREEN (existing impl).
"""

import gzip
import io
import subprocess
import sys
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import yaml

PROJECT_ROOT = Path(__file__).parent.parent


# =============================================================================
# 1. vaf_filter — _clean_bases (pure string processing)
# =============================================================================

class TestCleanBases:
    """
    _clean_bases strips samtools mpileup markup from a base string.
    Markup types: ^Q (read-start + mapq), $ (read-end),
                  -Nxxx (deletion), +Nxxx (insertion).
    """

    @pytest.fixture
    def clean_bases(self):
        from scripts.vaf_filter import _clean_bases
        return _clean_bases

    def test_plain_bases_unchanged(self, clean_bases):
        assert clean_bases("ACGT") == "ACGT"

    def test_ref_matches_unchanged(self, clean_bases):
        """Dots and commas (ref matches) are kept — callers count them."""
        assert clean_bases(".,.,") == ".,.,."[:-1]  # ".,." stripped of nothing
        assert "." in clean_bases(".,.ACGT")

    def test_strips_read_start_marker(self, clean_bases):
        """^Q marks read start followed by mapping quality char — both removed."""
        assert clean_bases("^!ACG") == "ACG"
        assert clean_bases("A^~CG") == "ACG"

    def test_strips_read_end_marker(self, clean_bases):
        """$ marks read end — removed."""
        assert clean_bases("ACG$T") == "ACGT"
        assert clean_bases("$ACG$") == "ACG"

    def test_strips_deletion(self, clean_bases):
        """-3ACG means deletion of 3 bases — the markup is removed."""
        result = clean_bases("A-3ACGt")
        assert "-" not in result
        assert result == "At"

    def test_strips_insertion(self, clean_bases):
        """+2AC means insertion of 2 bases — the markup is removed."""
        result = clean_bases("A+2ACt")
        assert "+" not in result
        assert result == "At"

    def test_strips_combined_markup(self, clean_bases):
        """Real mpileup string with multiple markup types."""
        raw = "^!ACG$T-2AG+3ACTt"
        result = clean_bases(raw)
        assert "^" not in result
        assert "$" not in result
        assert "-" not in result
        assert "+" not in result

    def test_empty_string(self, clean_bases):
        assert clean_bases("") == ""


# =============================================================================
# 2. vaf_filter — binom_pvalue (pure math)
# =============================================================================

class TestBinomPvalue:
    """binom_pvalue returns one-sided p-value for H0: VAF >= 0.5."""

    @pytest.fixture
    def binom_pvalue(self):
        from scripts.vaf_filter import binom_pvalue
        return binom_pvalue

    def test_zero_depth_returns_one(self, binom_pvalue):
        assert binom_pvalue(0, 0) == 1.0

    def test_somatic_mosaic_very_small_pvalue(self, binom_pvalue):
        """5/100 (5% VAF) is unambiguously somatic → p << 1e-6."""
        assert binom_pvalue(5, 100) < 1e-6

    def test_germline_het_large_pvalue(self, binom_pvalue):
        """50/100 (50% VAF) is germline het → p ≈ 0.5."""
        p = binom_pvalue(50, 100)
        assert p > 0.4

    def test_returns_float(self, binom_pvalue):
        p = binom_pvalue(10, 200)
        assert isinstance(p, float)

    def test_bounded_0_to_1(self, binom_pvalue):
        for alt, depth in [(0, 100), (50, 100), (100, 100), (5, 1000)]:
            p = binom_pvalue(alt, depth)
            assert 0.0 <= p <= 1.0, f"p={p} out of range for alt={alt}, depth={depth}"


# =============================================================================
# 3. vaf_filter — pileup_base_counts (mocked subprocess)
# =============================================================================

class TestPileupBaseCounts:
    """pileup_base_counts calls samtools mpileup via apptainer."""

    @pytest.fixture
    def pileup_base_counts(self):
        from scripts.vaf_filter import pileup_base_counts
        return pileup_base_counts

    def _mock_run(self, stdout: str):
        m = MagicMock()
        m.stdout = stdout
        return m

    def test_returns_depth_and_counts(self, pileup_base_counts):
        """Standard pileup output → correct depth and base counts."""
        mpileup_out = "chr1\t1000\tA\t4\tACGT\tIIII\n"
        with patch("subprocess.run", return_value=self._mock_run(mpileup_out)):
            depth, counts = pileup_base_counts("s.cram", "ref.fa", "chr1", "1000", "sam.sif")
        assert depth == 4
        assert counts["A"] == 1
        assert counts["C"] == 1

    def test_zero_coverage_returns_empty(self, pileup_base_counts):
        """Empty mpileup output (no coverage) → (0, {})."""
        with patch("subprocess.run", return_value=self._mock_run("")):
            depth, counts = pileup_base_counts("s.cram", "ref.fa", "chr1", "1000", "sam.sif")
        assert depth == 0
        assert counts == {}

    def test_subprocess_error_returns_empty(self, pileup_base_counts):
        """CalledProcessError → graceful (0, {}) with no exception."""
        with patch("subprocess.run", side_effect=subprocess.CalledProcessError(1, "apptainer")):
            depth, counts = pileup_base_counts("s.cram", "ref.fa", "chr1", "1000", "sam.sif")
        assert depth == 0
        assert counts == {}

    def test_counts_both_upper_and_lower(self, pileup_base_counts):
        """mpileup returns both upper (forward) and lower (reverse) strand bases."""
        mpileup_out = "chr1\t1000\tA\t4\tAaGg\tIIII\n"
        with patch("subprocess.run", return_value=self._mock_run(mpileup_out)):
            depth, counts = pileup_base_counts("s.cram", "ref.fa", "chr1", "1000", "sam.sif")
        assert counts["A"] == 1
        assert counts["a"] == 1
        assert counts["G"] == 1
        assert counts["g"] == 1


# =============================================================================
# 4. vaf_filter — filter_variants (mocked pileup)
# =============================================================================

class TestVafFilterVariants:
    """filter_variants integrates pileup_base_counts + passes_vaf_filter."""

    @pytest.fixture
    def filter_variants(self):
        from scripts.vaf_filter import filter_variants
        return filter_variants

    def test_somatic_variant_kept(self, filter_variants):
        """alt_n=10, depth=300 → somatic → kept in output."""
        lines = "chr1\t1000\tA\tG\n"
        infile = io.StringIO(lines)
        outfile = io.StringIO()
        with patch("scripts.vaf_filter.pileup_base_counts", return_value=(300, {"G": 10, "g": 0})):
            filter_variants(infile, "s.cram", "ref.fa", "sam.sif", outfile=outfile)
        assert "1000" in outfile.getvalue()

    def test_germline_variant_removed(self, filter_variants):
        """alt_n=50, depth=100 (VAF ~50%) → germline → removed."""
        lines = "chr1\t2000\tA\tT\n"
        infile = io.StringIO(lines)
        outfile = io.StringIO()
        with patch("scripts.vaf_filter.pileup_base_counts", return_value=(100, {"T": 50, "t": 0})):
            filter_variants(infile, "s.cram", "ref.fa", "sam.sif", outfile=outfile)
        assert "2000" not in outfile.getvalue()

    def test_zero_coverage_variant_removed(self, filter_variants):
        """No coverage (depth=0) → removed."""
        lines = "chr1\t3000\tA\tC\n"
        infile = io.StringIO(lines)
        outfile = io.StringIO()
        with patch("scripts.vaf_filter.pileup_base_counts", return_value=(0, {})):
            filter_variants(infile, "s.cram", "ref.fa", "sam.sif", outfile=outfile)
        assert "3000" not in outfile.getvalue()

    def test_comment_lines_skipped(self, filter_variants):
        """Lines starting with # are skipped (not passed through)."""
        lines = "#comment\nchr1\t1000\tA\tG\n"
        infile = io.StringIO(lines)
        outfile = io.StringIO()
        with patch("scripts.vaf_filter.pileup_base_counts", return_value=(300, {"G": 10})):
            filter_variants(infile, "s.cram", "ref.fa", "sam.sif", outfile=outfile)
        assert "#comment" not in outfile.getvalue()

    def test_short_lines_skipped(self, filter_variants):
        """Lines with fewer than 4 fields are silently skipped."""
        lines = "chr1\t1000\n"
        infile = io.StringIO(lines)
        outfile = io.StringIO()
        with patch("scripts.vaf_filter.pileup_base_counts", return_value=(100, {"A": 10})):
            filter_variants(infile, "s.cram", "ref.fa", "sam.sif", outfile=outfile)
        assert outfile.getvalue() == ""


# =============================================================================
# 5. vaf_filter — CLI (main)
# =============================================================================

class TestVafFilterCLI:
    """CLI contract tests for vaf_filter.py main()."""

    def test_help_exits_zero(self):
        result = subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "scripts" / "vaf_filter.py"), "--help"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert "--cram" in result.stdout

    def test_missing_required_args_exits_nonzero(self):
        """--cram is required; omitting it must exit non-zero."""
        result = subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "scripts" / "vaf_filter.py")],
            capture_output=True, text=True,
            input="",
        )
        assert result.returncode != 0

    def test_main_invocable_via_argparse(self, tmp_path):
        """main() runs when called with valid args (pileup mocked)."""
        from scripts.vaf_filter import main
        infile = tmp_path / "variants.txt"
        infile.write_text("chr1\t1000\tA\tG\n")
        argv = [
            "vaf_filter.py",
            str(infile),
            "--cram", "dummy.cram",
            "--ref", "dummy.fa",
            "--samtools-sif", "dummy.sif",
        ]
        with patch("sys.argv", argv), \
             patch("scripts.vaf_filter.pileup_base_counts", return_value=(300, {"G": 10})), \
             patch("sys.stdout", io.StringIO()):
            main()


# =============================================================================
# 6. pon_mask_filter — query_fasta (mocked subprocess)
# =============================================================================

class TestQueryFasta:
    """query_fasta calls samtools faidx via apptainer."""

    @pytest.fixture
    def query_fasta(self):
        from scripts.pon_mask_filter import query_fasta
        return query_fasta

    def test_returns_base_from_faidx_output(self, query_fasta):
        mock_result = MagicMock()
        mock_result.stdout = ">chr1:1000-1000\nA\n"
        with patch("subprocess.run", return_value=mock_result):
            base = query_fasta("pon.fa", "chr1", "1000", "sam.sif")
        assert base == "A"

    def test_returns_question_mark_on_error(self, query_fasta):
        with patch("subprocess.run", side_effect=subprocess.CalledProcessError(1, "apptainer")):
            base = query_fasta("pon.fa", "chr1", "1000", "sam.sif")
        assert base == "?"

    def test_returns_question_mark_on_empty_output(self, query_fasta):
        mock_result = MagicMock()
        mock_result.stdout = ">chr1:1000-1000\n"  # header only, no sequence
        with patch("subprocess.run", return_value=mock_result):
            base = query_fasta("pon.fa", "chr1", "1000", "sam.sif")
        assert base == "?"


# =============================================================================
# 7. pon_mask_filter — filter_variants (mocked query_fasta)
# =============================================================================

class TestPonMaskFilterVariants:
    """filter_variants integrates query_fasta + pon_passes."""

    @pytest.fixture
    def filter_variants(self):
        from scripts.pon_mask_filter import filter_variants
        return filter_variants

    def test_variant_not_in_pon_passes(self, filter_variants):
        """PON base = '*' (not in PON) → variant passes through."""
        lines = "chr1\t1000\tA\tG\n"
        infile = io.StringIO(lines)
        outfile = io.StringIO()
        with patch("scripts.pon_mask_filter.query_fasta", return_value="*"):
            filter_variants(infile, "pon.fa", "sam.sif", outfile=outfile)
        assert "1000" in outfile.getvalue()

    def test_exact_match_masked(self, filter_variants):
        """PON base = 'G', alt = 'G' → masked (removed)."""
        lines = "chr1\t2000\tA\tG\n"
        infile = io.StringIO(lines)
        outfile = io.StringIO()
        with patch("scripts.pon_mask_filter.query_fasta", return_value="G"):
            filter_variants(infile, "pon.fa", "sam.sif", outfile=outfile)
        assert "2000" not in outfile.getvalue()

    def test_iupac_match_masked(self, filter_variants):
        """PON base = 'R' (A|G), alt = 'G' → masked."""
        lines = "chr1\t3000\tA\tG\n"
        infile = io.StringIO(lines)
        outfile = io.StringIO()
        with patch("scripts.pon_mask_filter.query_fasta", return_value="R"):
            filter_variants(infile, "pon.fa", "sam.sif", outfile=outfile)
        assert "3000" not in outfile.getvalue()

    def test_iupac_non_match_passes(self, filter_variants):
        """PON base = 'R' (A|G), alt = 'C' → passes (not in R set)."""
        lines = "chr1\t4000\tA\tC\n"
        infile = io.StringIO(lines)
        outfile = io.StringIO()
        with patch("scripts.pon_mask_filter.query_fasta", return_value="R"):
            filter_variants(infile, "pon.fa", "sam.sif", outfile=outfile)
        assert "4000" in outfile.getvalue()

    def test_query_error_conservative_fail(self, filter_variants):
        """query_fasta returns '?' (error) → unknown code → conservative Fail → removed."""
        lines = "chr1\t5000\tA\tG\n"
        infile = io.StringIO(lines)
        outfile = io.StringIO()
        with patch("scripts.pon_mask_filter.query_fasta", return_value="?"):
            filter_variants(infile, "pon.fa", "sam.sif", outfile=outfile)
        assert "5000" not in outfile.getvalue()

    def test_comment_lines_skipped(self, filter_variants):
        lines = "#header\nchr1\t1000\tA\tG\n"
        infile = io.StringIO(lines)
        outfile = io.StringIO()
        with patch("scripts.pon_mask_filter.query_fasta", return_value="*"):
            filter_variants(infile, "pon.fa", "sam.sif", outfile=outfile)
        assert "#header" not in outfile.getvalue()


# =============================================================================
# 8. pon_mask_filter — CLI
# =============================================================================

class TestPonMaskFilterCLI:

    def test_help_exits_zero(self):
        result = subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "scripts" / "pon_mask_filter.py"), "--help"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert "--pon-fasta" in result.stdout

    def test_missing_required_args_exits_nonzero(self):
        result = subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "scripts" / "pon_mask_filter.py")],
            capture_output=True, text=True, input="",
        )
        assert result.returncode != 0

    def test_main_invocable_via_argparse(self, tmp_path):
        from scripts.pon_mask_filter import main
        infile = tmp_path / "variants.txt"
        infile.write_text("chr1\t1000\tA\tG\n")
        argv = [
            "pon_mask_filter.py",
            str(infile),
            "--pon-fasta", "dummy.fa",
            "--samtools-sif", "dummy.sif",
        ]
        with patch("sys.argv", argv), \
             patch("scripts.pon_mask_filter.query_fasta", return_value="*"), \
             patch("sys.stdout", io.StringIO()):
            main()


# =============================================================================
# 9. germline_filter — CLI (main)
# =============================================================================

class TestGermlineFilterCLI:

    def test_help_exits_zero(self):
        result = subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "scripts" / "germline_filter.py"), "--help"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert "--variants" in result.stdout

    def test_missing_variants_flag_exits_nonzero(self):
        """--variants is required."""
        result = subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "scripts" / "germline_filter.py")],
            capture_output=True, text=True, input="chr1\t100\tA\tG\n",
        )
        assert result.returncode != 0

    def test_filters_via_stdin(self, tmp_path):
        """End-to-end: germline variant on stdin is removed."""
        gz = tmp_path / "gnomad.txt.gz"
        with gzip.open(gz, "wt") as fh:
            fh.write("1\t100\tA\tG\n")

        result = subprocess.run(
            [
                sys.executable,
                str(PROJECT_ROOT / "scripts" / "germline_filter.py"),
                "--variants", str(gz),
            ],
            input="chr1\t100\tA\tG\nchr1\t200\tC\tT\n",
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert "100" not in result.stdout   # germline removed
        assert "200" in result.stdout       # somatic kept

    def test_main_invocable_via_argparse(self, tmp_path):
        from scripts.germline_filter import main
        gz = tmp_path / "gnomad.txt.gz"
        with gzip.open(gz, "wt") as fh:
            fh.write("1\t100\tA\tG\n")
        infile = tmp_path / "variants.txt"
        infile.write_text("chr1\t100\tA\tG\nchr1\t200\tC\tT\n")
        argv = ["germline_filter.py", str(infile), "--variants", str(gz)]
        with patch("sys.argv", argv), patch("sys.stdout", io.StringIO()):
            main()


# =============================================================================
# 10. prepare_containers — function and CLI
# =============================================================================

class TestPrepareContainers:

    @pytest.fixture
    def containers_yaml(self, tmp_path):
        """Minimal containers.yaml with two tools, one SIF existing."""
        sif_exists = tmp_path / "containers" / "existing.sif"
        sif_exists.parent.mkdir()
        sif_exists.touch()

        cfg = {
            "tool_a": {
                "name": "tool_a", "version": "1.0",
                "uri": "docker://example/tool_a:1.0",
                "sif": str(sif_exists),
            },
            "tool_b": {
                "name": "tool_b", "version": "2.0",
                "uri": "docker://example/tool_b:2.0",
                "sif": str(tmp_path / "containers" / "missing.sif"),
            },
        }
        yaml_path = tmp_path / "containers.yaml"
        with open(yaml_path, "w") as fh:
            yaml.dump(cfg, fh)
        return str(yaml_path)

    def test_check_only_does_not_pull(self, containers_yaml, tmp_path, capsys):
        """--check-only prints warning, never calls apptainer pull."""
        from scripts.prepare_containers import prepare_containers
        with patch("subprocess.run") as mock_run:
            prepare_containers(containers_yaml, pull=False)
            mock_run.assert_not_called()
        captured = capsys.readouterr()
        assert "WARNING" in captured.out or "missing.sif" in captured.out

    def test_existing_sif_skipped(self, containers_yaml, tmp_path):
        """SIF that already exists is not pulled again."""
        from scripts.prepare_containers import prepare_containers
        pulled = []
        def fake_run(cmd, **kw):
            pulled.append(cmd)
        with patch("subprocess.run", side_effect=fake_run):
            prepare_containers(containers_yaml, pull=True)
        # Only tool_b (missing) should be pulled
        assert len(pulled) == 1
        assert "tool_b" in " ".join(pulled[0]) or "missing" in " ".join(pulled[0])

    def test_all_present_prints_message(self, tmp_path, capsys):
        """When all SIFs exist, prints 'All SIF files present.'"""
        from scripts.prepare_containers import prepare_containers
        sif = tmp_path / "containers" / "all.sif"
        sif.parent.mkdir()
        sif.touch()
        cfg = {"tool": {"name": "t", "version": "1", "uri": "docker://x:1", "sif": str(sif)}}
        yaml_path = tmp_path / "all.yaml"
        with open(yaml_path, "w") as fh:
            yaml.dump(cfg, fh)
        prepare_containers(str(yaml_path), pull=False)
        assert "All SIF files present" in capsys.readouterr().out


class TestPrepareContainersCLI:

    def test_help_exits_zero(self):
        result = subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "scripts" / "prepare_containers.py"), "--help"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0

    def test_check_only_flag_does_not_pull(self, tmp_path):
        """--check-only must not invoke apptainer."""
        sif = tmp_path / "containers" / "s.sif"
        sif.parent.mkdir()
        sif.touch()
        cfg = {"t": {"name": "t", "version": "1", "uri": "docker://x:1", "sif": str(sif)}}
        yaml_path = tmp_path / "c.yaml"
        with open(yaml_path, "w") as fh:
            yaml.dump(cfg, fh)

        result = subprocess.run(
            [
                sys.executable,
                str(PROJECT_ROOT / "scripts" / "prepare_containers.py"),
                str(yaml_path),
                "--check-only",
            ],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert "All SIF files present" in result.stdout


# =============================================================================
# 11. auto_params — CLI (main)
# =============================================================================

class TestAutoParamsCLI:

    def test_help_exits_zero(self):
        result = subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "scripts" / "auto_params.py"), "--help"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0
        assert "--output" in result.stdout

    def test_missing_input_exits_nonzero(self):
        result = subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "scripts" / "auto_params.py")],
            capture_output=True, text=True,
        )
        assert result.returncode != 0

    def test_cli_writes_resolved_params(self, tmp_path):
        """CLI with a real FASTQ produces a valid resolved_params.yaml."""
        # Reuse _make_fastq helper via inline creation
        fq = tmp_path / "novaseq.R1.fastq.gz"
        import gzip as _gz
        with _gz.open(fq, "wt") as fh:
            for i in range(5):
                fh.write(f"@A00100:123:AABBCCDD:1:1101:{i}:{i} 1:N:0:\n")
                fh.write("ACGT\n+\nIIII\n")
        out_yaml = tmp_path / "resolved.yaml"

        result = subprocess.run(
            [
                sys.executable,
                str(PROJECT_ROOT / "scripts" / "auto_params.py"),
                str(fq),
                "--output", str(out_yaml),
            ],
            capture_output=True, text=True,
            cwd=str(PROJECT_ROOT),
        )
        assert result.returncode == 0, f"STDERR: {result.stderr}"
        assert out_yaml.exists()
        with open(out_yaml) as fh:
            data = yaml.safe_load(fh)
        assert data["sequencer"] == "NovaSeq 6000"
        assert data["optical_duplicate_pixel_distance"] == 2500
        assert "sequencer_evidence" in data

    def test_cli_prints_summary(self, tmp_path):
        """CLI prints sequencer, ODPD, and output path to stdout."""
        fq = tmp_path / "test.R1.fastq.gz"
        import gzip as _gz
        with _gz.open(fq, "wt") as fh:
            fh.write("@A00100:1:FC:1:1:1:1 1:N:0:\nACGT\n+\nIIII\n")
        out_yaml = tmp_path / "r.yaml"

        result = subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "scripts" / "auto_params.py"),
             str(fq), "--output", str(out_yaml)],
            capture_output=True, text=True, cwd=str(PROJECT_ROOT),
        )
        assert "NovaSeq 6000" in result.stdout
        assert "2500" in result.stdout
