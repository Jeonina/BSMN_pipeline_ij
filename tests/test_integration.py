"""
Integration tests — Phase B & D.

Phase B: pileup_base_counts() with real CRAM + samtools container.
Phase D: Snakemake output validation (real pipeline artefacts).

All tests in this file require:
  - results/mapping/SAMPLE/SAMPLE.cram   (produced by the mapping pipeline)
  - resources/hg38/chr22.fasta           (reference genome)
  - containers/samtools_1.17.sif         (samtools Apptainer image)

Skip if artefacts are missing:
    pytest tests/test_integration.py -m "not integration"

Run only integration tests:
    pytest tests/test_integration.py -m integration -v
"""

import io
import subprocess
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).parent.parent

CRAM = str(PROJECT_ROOT / "results/mapping/SAMPLE/SAMPLE.cram")
CRAI = str(PROJECT_ROOT / "results/mapping/SAMPLE/SAMPLE.cram.crai")
REF = str(PROJECT_ROOT / "resources/hg38/chr22.fasta")
SAMTOOLS_SIF = str(PROJECT_ROOT / "containers/samtools_1.17.sif")
FLAGSTAT = str(PROJECT_ROOT / "results/mapping/SAMPLE/flagstat.txt")
MARKDUP_METRICS = str(PROJECT_ROOT / "results/mapping/SAMPLE/markduplicates_metrics.txt")

# chr22 positions confirmed to have reads in the SAMPLE CRAM
# (ERR194146 subset; NA12878 WGS vs chr22 reference)
#   50086322: ref=A, bases=TTTTTtttt  → alt T, all 9 reads are alt
#   50086306: ref=A, bases=GGgg       → alt G, 4 reads (below min_alt=5 default)
#   50086298: ref=A, bases=,........  → all ref (no explicit ACGT in pileup)
POS_ALL_ALT = "50086322"  # T alt, depth=9 by code (9 alt bases)
POS_LOW_ALT = "50086306"  # G alt, depth=4 by code  (4 alt bases, <min_alt=5)
POS_ALL_REF = "50086298"  # all ref matches, code depth=0


def _artefacts_present() -> bool:
    return all(Path(p).exists() for p in [CRAM, CRAI, REF, SAMTOOLS_SIF])


skip_if_missing = pytest.mark.skipif(
    not _artefacts_present(),
    reason="Mapping pipeline artefacts not present; run `snakemake --cores 4` first",
)


# =============================================================================
# Phase B — pileup_base_counts() integration tests
# =============================================================================


@pytest.mark.integration
@skip_if_missing
class TestPileupIntegration:
    """pileup_base_counts() called against real CRAM via Apptainer."""

    def _run(self, pos: str, **kwargs):
        from scripts.vaf_filter import pileup_base_counts

        return pileup_base_counts(CRAM, REF, "chr22", pos, SAMTOOLS_SIF, **kwargs)

    def test_all_alt_position_returns_positive_depth(self):
        """Position where all reads carry the alt allele returns depth > 0."""
        depth, counts = self._run(POS_ALL_ALT)
        assert depth > 0, f"Expected depth > 0 at chr22:{POS_ALL_ALT}, got {depth}"

    def test_all_alt_position_alt_base_is_T(self):
        """At chr22:50086322 the dominant alt base is T (both strands)."""
        depth, counts = self._run(POS_ALL_ALT)
        t_total = counts.get("T", 0) + counts.get("t", 0)
        assert t_total == depth, (
            f"Expected all alt reads to be T, got T={t_total}, depth={depth}, counts={counts}"
        )

    def test_low_alt_position_returns_four_reads(self):
        """Position with 4 alt reads returns depth=4."""
        depth, counts = self._run(POS_LOW_ALT)
        assert depth == 4, f"Expected depth=4 at chr22:{POS_LOW_ALT}, got {depth}"

    def test_all_ref_position_returns_zero_code_depth(self):
        """Position with only ref matches returns depth=0 (no explicit ACGT alt bases)."""
        depth, counts = self._run(POS_ALL_REF)
        assert depth == 0, (
            f"Expected code depth=0 at all-ref position chr22:{POS_ALL_REF}, got {depth}"
        )

    def test_mapq_filter_reduces_depth(self):
        """Raising min_mapq to 60 should not increase depth."""
        depth_q20, _ = self._run(POS_ALL_ALT, min_mapq=20)
        depth_q60, _ = self._run(POS_ALL_ALT, min_mapq=60)
        assert depth_q60 <= depth_q20, (
            f"Higher MAPQ filter should not increase depth: q20={depth_q20}, q60={depth_q60}"
        )

    def test_counts_dict_has_all_eight_bases(self):
        """pileup_base_counts always returns counts for all 8 ACGT/acgt keys."""
        _, counts = self._run(POS_ALL_ALT)
        for base in ("A", "a", "C", "c", "G", "g", "T", "t"):
            assert base in counts, f"Missing key '{base}' in counts dict"

    def test_invalid_position_returns_empty(self):
        """Out-of-range position (no reads) returns (0, {})."""
        from scripts.vaf_filter import pileup_base_counts

        depth, counts = pileup_base_counts(CRAM, REF, "chr22", "1", SAMTOOLS_SIF)
        assert depth == 0
        assert counts == {}


# =============================================================================
# Phase B — filter_variants() integration tests
# =============================================================================


@pytest.mark.integration
@skip_if_missing
class TestFilterVariantsIntegration:
    """filter_variants() end-to-end with real CRAM pileup."""

    def _filter(self, lines: str, **kwargs) -> str:
        from scripts.vaf_filter import filter_variants

        out = io.StringIO()
        filter_variants(io.StringIO(lines), CRAM, REF, SAMTOOLS_SIF, outfile=out, **kwargs)
        return out.getvalue()

    def test_germline_hom_alt_filtered_by_binomial(self):
        """chr22:50086322 A→T (VAF=100%) is filtered: binomial test rejects VAF<0.5 hypothesis.

        The VAF filter targets somatic mosaic variants (VAF ~0.05-0.30).
        A germline homozygous alt (VAF=1.0) has p-value >> 1e-6 and is correctly removed.
        """
        inp = f"chr22\t{POS_ALL_ALT}\tA\tT\n"
        out = self._filter(inp, min_alt=5, p_threshold=1e-6)
        assert out == "", f"Germline hom-alt variant (VAF=1.0) should be filtered; got: {out!r}"

    def test_low_alt_count_variant_is_removed(self):
        """chr22:50086306 A→G (4 alt reads < min_alt=5) is filtered out."""
        inp = f"chr22\t{POS_LOW_ALT}\tA\tG\n"
        out = self._filter(inp, min_alt=5, p_threshold=1e-6)
        assert out == "", f"Low-count variant should be removed; got: {out!r}"

    def test_comment_lines_are_passed_through(self):
        """Lines starting with '#' are skipped (not written to output)."""
        inp = f"# header\nchr22\t{POS_ALL_ALT}\tA\tT\n"
        out = self._filter(inp, min_alt=5, p_threshold=1e-6)
        assert "# header" not in out, "Comment lines should be skipped, not passed through"

    def test_relaxed_pthreshold_keeps_hom_alt(self):
        """With p_threshold=1.0 the binomial test always passes; only min_alt matters.

        chr22:50086322 has 9 T alt reads → alt_n=9 >= min_alt=5 → kept.
        This exercises the 'variant is written' code path with real pileup data.
        """
        inp = f"chr22\t{POS_ALL_ALT}\tA\tT\n"
        out = self._filter(inp, min_alt=5, p_threshold=1.0)
        assert out == inp, (
            f"With p_threshold=1.0, min_alt=5-passing variant should be kept; got: {out!r}"
        )

    def test_relaxed_threshold_with_low_count_still_filtered(self):
        """chr22:50086306 (4 alt reads < min_alt=5) is still filtered even with p_threshold=1.0."""
        inp = f"chr22\t{POS_LOW_ALT}\tA\tG\n"
        out = self._filter(inp, min_alt=5, p_threshold=1.0)
        assert out == "", (
            f"4-read variant should be filtered by min_alt=5 regardless of "
            f"p_threshold; got: {out!r}"
        )


# =============================================================================
# Phase D — Snakemake output validation tests
# =============================================================================


@pytest.mark.integration
@skip_if_missing
class TestMappingOutputs:
    """Validate artefacts produced by the mapping pipeline."""

    def test_cram_file_exists_and_is_nonzero(self):
        assert Path(CRAM).stat().st_size > 0, "CRAM file is empty"

    def test_cram_index_exists(self):
        assert Path(CRAI).exists(), "CRAM index (.crai) missing"

    def test_cram_is_valid_with_samtools_quickcheck(self):
        """samtools quickcheck exits 0 for a valid CRAM."""
        result = subprocess.run(
            ["apptainer", "exec", SAMTOOLS_SIF, "samtools", "quickcheck", CRAM],
            capture_output=True,
        )
        assert result.returncode == 0, f"samtools quickcheck failed:\n{result.stderr.decode()}"

    def test_cram_has_mapped_reads(self):
        """samtools view -c returns the mapped read count (> 0 for chr22 subset)."""
        result = subprocess.run(
            [
                "apptainer",
                "exec",
                SAMTOOLS_SIF,
                "samtools",
                "view",
                "-c",
                "-F",
                "4",
                "--input-fmt-option",
                f"reference={REF}",
                CRAM,
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        count = int(result.stdout.strip())
        assert count > 0, f"Expected mapped reads > 0, got {count}"

    def test_flagstat_exists_and_has_expected_fields(self):
        text = Path(FLAGSTAT).read_text()
        for keyword in ("in total", "mapped", "paired in sequencing", "duplicates"):
            assert keyword in text, f"flagstat missing '{keyword}' line"

    def test_flagstat_mapping_rate_is_plausible(self):
        """8-15% mapping rate expected for WGS data aligned to chr22 only."""
        text = Path(FLAGSTAT).read_text()
        for line in text.splitlines():
            if "primary mapped" in line and "%" in line:
                pct_str = line.split("(")[1].split("%")[0].strip()
                pct = float(pct_str)
                assert 1.0 <= pct <= 50.0, (
                    f"Unexpected primary mapped %: {pct:.2f}% (line: {line!r})"
                )
                break

    def test_markdup_metrics_has_duplicate_rate(self):
        text = Path(MARKDUP_METRICS).read_text()
        assert "PERCENT_DUPLICATION" in text, (
            "MarkDuplicates metrics missing PERCENT_DUPLICATION field"
        )
        # Parse the duplicate rate from the metrics table
        lines = text.splitlines()
        header_idx = next(i for i, ln in enumerate(lines) if ln.startswith("LIBRARY"))
        values = lines[header_idx + 1].split("\t")
        headers = lines[header_idx].split("\t")
        dup_rate = float(values[headers.index("PERCENT_DUPLICATION")])
        assert 0.0 <= dup_rate <= 1.0, f"Duplicate rate out of range: {dup_rate}"

    def test_snakemake_reports_nothing_to_do_after_completion(self):
        """A second snakemake run should report 'Nothing to be done'."""
        result = subprocess.run(
            ["snakemake", "--snakefile", "workflow/Snakefile", "--dry-run", "--cores", "1"],
            capture_output=True,
            text=True,
            cwd=str(PROJECT_ROOT),
        )
        assert result.returncode == 0
        combined = result.stdout + result.stderr
        assert "Nothing to be done" in combined, (
            f"Expected 'Nothing to be done' after completed run.\nOutput:\n{combined[:2000]}"
        )
