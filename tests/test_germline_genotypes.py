"""Tests for scripts/germline_genotypes.py (no bcftools needed)."""

import pytest

from scripts.germline_genotypes import af_to_dosage, cmd_merge, merge_genotypes, parse_query


def test_af_bands():
    bands = (0.25, 0.75, 0.85)
    assert [af_to_dosage(x, *bands) for x in (0.03, 0.24, 0.25, 0.5, 0.75, 0.8, 0.85, 1.0)] == \
        [None, None, "1", "1", "1", ".", "2", "2"]


def test_parse_query_keeps_biallelic_snps_on_main_chroms():
    text = ("chr1\t10\tA\tG\t0.48\n"          # het
            "chr1\t20\tA\tG,T\t0.4,0.5\n"     # multiallelic: skipped
            "chr1\t30\tAT\tA\t0.5\n"          # indel: skipped
            "chrUn_x\t5\tA\tC\t0.5\n"         # off main chromosomes: skipped
            "chrX\t40\tC\tT\t0.99\n"          # hom-alt
            "chr2\t50\tC\tT\t0.04\n"          # noise / contamination: not carried
            "chr2\t60\tC\tT\t.\n")
    assert parse_query(text, 0.25, 0.75, 0.85) == [(("chr1", 10, "A", "G"), "1"), (("chrX", 40, "C", "T"), "2")]


def test_merge_fills_hom_ref_and_drops_rare_sites():
    s1 = [(("chr1", 9, "C", "T"), "2"), (("chr2", 5, "A", "G"), "1")]
    s2 = [(("chr1", 9, "C", "T"), "."), (("chr2", 5, "A", "G"), "2")]
    s3 = [(("chr1", 3, "G", "A"), "1")]                          # seen once: dropped at min 2
    out = list(merge_genotypes([s1, s2, s3], min_carriers=2))
    assert out == [(("chr1", 9, "C", "T"), ["1/1", "./.", "0/0"]),
                   (("chr2", 5, "A", "G"), ["0/1", "1/1", "0/0"])]


def test_merge_cli_writes_vcf_and_sex_file(tmp_path):
    info = tmp_path / "info.tsv"
    info.write_text("sample_id\tbr_region\tsex\nS1\tBrA_DLPFC\tM\nS2\tBrB_Caudate\tF\nS3\tBrC_DLPFC\tM\n")
    (tmp_path / "S1.gt.tsv").write_text("chr1\t9\tC\tT\t1\n")
    (tmp_path / "S2.gt.tsv").write_text("chr1\t9\tC\tT\t2\n")
    (tmp_path / "S3.gt.tsv").write_text("chr1\t9\tC\tT\t1\n")   # excluded below

    class A:
        sample_info, out_dir, exclude, min_carriers = str(info), str(tmp_path), "BrC_DLPFC", 1

    assert cmd_merge(A) == 0
    lines = (tmp_path / "cohort_germline.vcf").read_text().splitlines()
    assert lines[-2].endswith("BrA_DLPFC\tBrB_Caudate")
    assert lines[-1] == "chr1\t9\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\t1/1"
    assert (tmp_path / "sex_update.txt").read_text() == "BrA_DLPFC\tBrA_DLPFC\t1\nBrB_Caudate\tBrB_Caudate\t2\n"


def test_merge_keeps_alleles_at_one_position_apart_and_orders_chrx_last():
    s1 = [(("chr2", 7, "A", "T"), "1"), (("chr2", 7, "A", "G"), "1"), (("chrX", 1, "C", "T"), "2")]
    s2 = [(("chr2", 7, "A", "G"), "2"), (("chrX", 1, "C", "T"), "1")]
    out = list(merge_genotypes([s1, s2], min_carriers=1))
    assert out == [(("chr2", 7, "A", "G"), ["0/1", "1/1"]),
                   (("chr2", 7, "A", "T"), ["0/1", "0/0"]),
                   (("chrX", 1, "C", "T"), ["1/1", "0/1"])]


def test_merge_refuses_unsorted_input():
    s1 = [(("chr2", 5, "A", "G"), "1"), (("chr1", 9, "C", "T"), "2")]
    with pytest.raises(ValueError, match="genome order"):
        list(merge_genotypes([s1], min_carriers=1))
