"""Tests for scripts/contamination_source.py (no bcftools needed)."""

from scripts.contamination_source import cmd_analyze, load_final, parse_query, source_table


def test_parse_query_keeps_germline_alleles_only():
    text = "chr1\t100\tG\tA\t0.48\nchr1\t200\tT\tC,G\t0.02,0.51\nchr2\t5\tA\tC\t.\nchr2\t9\tA\tC\t0.05\n"
    assert parse_query(text, 0.3) == {("chr1", 100, "G", "A"), ("chr1", 200, "T", "G")}


def test_source_table_finds_foreign_donor_and_skips_own_donor():
    donor_of = {"S1": "BrA", "S2": "BrA", "S3": "BrB", "S4": "BrC"}
    calls = {("chr1", i, "A", "G") for i in range(10)}
    finals = {"S1": calls, "S2": set(), "S3": set(), "S4": set()}
    germ = {
        "S1": set(),
        "S2": set(calls),                                   # own donor: shared calls are real, not contamination
        "S3": {("chr1", i, "A", "G") for i in range(8)},    # BrB carries 8 of S1's calls as germline
        "S4": {("chr1", 9, "A", "G"), ("chr1", 0, "A", "T")},  # wrong allele at pos 0 does not count
    }
    row = next(r for r in source_table(finals, germ, donor_of) if r["sample"] == "S1")
    assert (row["top_source"], row["top_hits"], row["top_frac"]) == ("BrB", 8, 0.8)
    assert (row["second_source"], row["second_hits"], row["any_foreign_germline"]) == ("BrC", 1, 9)
    empty = next(r for r in source_table(finals, germ, donor_of) if r["sample"] == "S3")
    assert (empty["top_source"], empty["top_frac"]) == ("", 0.0)


def test_analyze_writes_table(tmp_path):
    info = tmp_path / "info.tsv"
    info.write_text("sample_id\tbr_region\tsex\nS1\tBrA_DLPFC\tM\nS3\tBrB_Caudate\tF\n")
    for s, lines in {"S1": "chr1\t1\tA\tG\nchr1\t2\tC\tT\n", "S3": ""}.items():
        (tmp_path / "filtering" / s).mkdir(parents=True)
        (tmp_path / "filtering" / s / f"{s}.final.txt").write_text(lines)
    out = tmp_path / "out"
    out.mkdir()
    (out / "S1.germline.tsv").write_text("")
    (out / "S3.germline.tsv").write_text("chr1\t1\tA\tG\n")

    class A:
        sample_info, results, out_dir = str(info), str(tmp_path), str(out)

    assert cmd_analyze(A) == 0
    rows = (out / "contamination_source.tsv").read_text().splitlines()
    assert rows[1].split("\t")[:7] == ["S1", "BrA_DLPFC", "BrA", "2", "BrB", "1", "0.5"]
    assert load_final(str(out / "S3.germline.tsv")) == {("chr1", 1, "A", "G")}
