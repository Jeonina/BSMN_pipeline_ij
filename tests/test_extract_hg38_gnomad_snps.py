"""Tests for scripts/extract_hg38_gnomad_snps.py (M-FIX-004).

The script extracts SNPs above an AF threshold from the hg38 af-only gnomAD
VCF (Mutect2 germline_resource) into the lookup table format consumed by
scripts/germline_filter.py.

Contract:
  Input:  gzipped VCF (hg38 coordinates, possibly multiallelic)
  Output: gzipped TSV — bare-chrom \t pos \t ref \t alt (one per allele)
  Skip:   indels, multi-char alts, AF <= threshold, malformed AF

RED-GREEN-REFACTOR — tests written before implementation.
"""

from __future__ import annotations

import gzip
import subprocess
import sys
import time
from collections.abc import Iterable
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).parent.parent
SCRIPT = PROJECT_ROOT / "scripts" / "extract_hg38_gnomad_snps.py"


def write_vcf_gz(path: Path, body_lines: Iterable[str]) -> None:
    """Write a minimal gzipped VCF with a fixed header plus body_lines."""
    header = [
        "##fileformat=VCFv4.2",
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    with gzip.open(path, "wt") as fh:
        for line in header:
            fh.write(line + "\n")
        for line in body_lines:
            fh.write(line + "\n")


def read_output(path: Path) -> list[tuple[str, ...]]:
    """Read gzipped TSV output as a list of tuple-rows."""
    rows: list[tuple[str, ...]] = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            rows.append(tuple(line.split("\t")))
    return rows


def run_cli(*args: str) -> subprocess.CompletedProcess[str]:
    """Invoke the CLI as a subprocess."""
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True,
        text=True,
        cwd=str(PROJECT_ROOT),
    )


# =============================================================================
# Unit-level tests against the importable functions
# =============================================================================


class TestSingleAlleleSelection:
    """Filtering decisions for biallelic records."""

    @pytest.fixture
    def parse(self):
        from scripts.extract_hg38_gnomad_snps import iter_passing_alleles

        return iter_passing_alleles

    def test_biallelic_snv_above_threshold_emitted(self, parse):
        line = "chr1\t100\t.\tA\tC\t.\t.\tAF=0.05"
        result = list(parse(line, af_threshold=0.001, strip_chr=True))
        assert result == [("1", "100", "A", "C")]

    def test_biallelic_snv_at_or_below_threshold_skipped(self, parse):
        line = "chr1\t100\t.\tA\tC\t.\t.\tAF=0.001"
        assert list(parse(line, af_threshold=0.001, strip_chr=True)) == []
        line2 = "chr1\t100\t.\tA\tC\t.\t.\tAF=0.0001"
        assert list(parse(line2, af_threshold=0.001, strip_chr=True)) == []

    def test_indel_ref_skipped(self, parse):
        line = "chr1\t100\t.\tAT\tA\t.\t.\tAF=0.5"
        assert list(parse(line, af_threshold=0.001, strip_chr=True)) == []

    def test_multi_char_alt_skipped(self, parse):
        line = "chr1\t300\t.\tA\tCC\t.\t.\tAF=0.5"
        assert list(parse(line, af_threshold=0.001, strip_chr=True)) == []

    def test_no_strip_chr_keeps_prefix(self, parse):
        line = "chr1\t100\t.\tA\tC\t.\t.\tAF=0.05"
        result = list(parse(line, af_threshold=0.001, strip_chr=False))
        assert result == [("chr1", "100", "A", "C")]


class TestMultiallelic:
    """Multiallelic records: per-allele filtering with index-matched AF."""

    @pytest.fixture
    def parse(self):
        from scripts.extract_hg38_gnomad_snps import iter_passing_alleles

        return iter_passing_alleles

    def test_three_alts_mixed_af(self, parse):
        # A>C passes (0.005), A>G fails (0.0005), A>T passes (0.5)
        line = "chr1\t100\t.\tA\tC,G,T\t.\t.\tAF=0.005,0.0005,0.5"
        result = list(parse(line, af_threshold=0.001, strip_chr=True))
        assert result == [
            ("1", "100", "A", "C"),
            ("1", "100", "A", "T"),
        ]

    def test_multiallelic_mixed_with_indel_alt(self, parse):
        # 2nd alt is an indel (CC), must be skipped; SNVs still considered
        line = "chr1\t100\t.\tA\tC,CC,T\t.\t.\tAF=0.5,0.5,0.5"
        result = list(parse(line, af_threshold=0.001, strip_chr=True))
        assert result == [
            ("1", "100", "A", "C"),
            ("1", "100", "A", "T"),
        ]


class TestMalformedInput:
    """Malformed records must be skipped gracefully, never crash."""

    @pytest.fixture
    def parse(self):
        from scripts.extract_hg38_gnomad_snps import iter_passing_alleles

        return iter_passing_alleles

    def test_missing_af_field_skipped(self, parse):
        line = "chr1\t100\t.\tA\tC\t.\t.\tDP=42"
        assert list(parse(line, af_threshold=0.001, strip_chr=True)) == []

    def test_malformed_af_value_skipped(self, parse):
        line = "chr1\t100\t.\tA\tC\t.\t.\tAF=."
        assert list(parse(line, af_threshold=0.001, strip_chr=True)) == []
        line2 = "chr1\t100\t.\tA\tC\t.\t.\tAF=notanumber"
        assert list(parse(line2, af_threshold=0.001, strip_chr=True)) == []

    def test_af_count_mismatch_with_alts_skipped(self, parse):
        # 3 alts but only 2 AF values — entire record is suspect, skip
        line = "chr1\t100\t.\tA\tC,G,T\t.\t.\tAF=0.5,0.5"
        assert list(parse(line, af_threshold=0.001, strip_chr=True)) == []


# =============================================================================
# CLI integration tests
# =============================================================================


class TestCli:
    def test_header_only_input_produces_empty_gz(self, tmp_path: Path):
        vcf = tmp_path / "in.vcf.gz"
        out = tmp_path / "out.txt.gz"
        write_vcf_gz(vcf, [])
        result = run_cli("--input", str(vcf), "--output", str(out), "--af-threshold", "0.001")
        assert result.returncode == 0, result.stderr
        assert out.exists()
        assert read_output(out) == []

    def test_mixed_records_emits_expected_rows(self, tmp_path: Path):
        vcf = tmp_path / "in.vcf.gz"
        out = tmp_path / "out.txt.gz"
        write_vcf_gz(
            vcf,
            [
                # multiallelic — A>C passes, A>G fails, A>T passes
                "chr1\t100\t.\tA\tC,G,T\t.\t.\tAF=0.005,0.0005,0.5",
                # biallelic SNV passes
                "chr1\t200\t.\tA\tC\t.\t.\tAF=0.05",
                # alt is indel — skip
                "chr1\t300\t.\tA\tCC\t.\t.\tAF=0.05",
                # ref is indel — skip
                "chr1\t400\t.\tAT\tA\t.\t.\tAF=0.5",
                # below threshold — skip
                "chr1\t500\t.\tG\tA\t.\t.\tAF=0.0005",
            ],
        )
        result = run_cli("--input", str(vcf), "--output", str(out), "--af-threshold", "0.001")
        assert result.returncode == 0, result.stderr
        rows = read_output(out)
        assert rows == [
            ("1", "100", "A", "C"),
            ("1", "100", "A", "T"),
            ("1", "200", "A", "C"),
        ]

    def test_no_strip_chr_flag_keeps_prefix(self, tmp_path: Path):
        vcf = tmp_path / "in.vcf.gz"
        out = tmp_path / "out.txt.gz"
        write_vcf_gz(vcf, ["chr1\t100\t.\tA\tC\t.\t.\tAF=0.05"])
        result = run_cli(
            "--input",
            str(vcf),
            "--output",
            str(out),
            "--af-threshold",
            "0.001",
            "--no-strip-chr",
        )
        assert result.returncode == 0, result.stderr
        assert read_output(out) == [("chr1", "100", "A", "C")]

    def test_missing_input_file_nonzero_exit(self, tmp_path: Path):
        out = tmp_path / "out.txt.gz"
        result = run_cli(
            "--input",
            str(tmp_path / "does_not_exist.vcf.gz"),
            "--output",
            str(out),
            "--af-threshold",
            "0.001",
        )
        assert result.returncode != 0
        assert not out.exists() or out.stat().st_size == 0

    def test_speed_regression_under_one_second(self, tmp_path: Path):
        vcf = tmp_path / "in.vcf.gz"
        out = tmp_path / "out.txt.gz"
        body = [f"chr1\t{1000 + i}\t.\tA\tC\t.\t.\tAF=0.05" for i in range(1000)]
        write_vcf_gz(vcf, body)
        t0 = time.perf_counter()
        result = run_cli("--input", str(vcf), "--output", str(out), "--af-threshold", "0.001")
        elapsed = time.perf_counter() - t0
        assert result.returncode == 0, result.stderr
        assert len(read_output(out)) == 1000
        # Process startup dominates; the parse itself is ~milliseconds.
        assert elapsed < 1.0, f"Processing 1000 records took {elapsed:.3f}s"
