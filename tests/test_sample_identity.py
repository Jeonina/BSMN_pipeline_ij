"""Tests for scripts/sample_identity.py (no bcftools/samtools needed)."""

import csv
import gzip
import random

from scripts.sample_identity import infer_sex, main, parse_mpileup, write_sites


def test_write_sites_keeps_common_biallelic_snps_outside_par(tmp_path):
    vcf = tmp_path / "s.vcf.gz"
    with gzip.open(vcf, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, ref, alt, af in [
            ("chr22", 100, "A", "G", 0.3),      # kept
            ("chr22", 200, "A", "G", 0.01),     # too rare
            ("chr22", 300, "AT", "A", 0.3),     # indel
            ("chr21", 400, "A", "G", 0.3),      # other chromosome
            ("chrX", 500_000, "C", "T", 0.4),   # PAR1
            ("chrX", 5_000_000, "C", "T", 0.4), # kept
        ]:
            fh.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tAF={af}\n")
    out = tmp_path / "sites.tsv"
    assert write_sites(str(vcf), str(out)) == 2
    assert out.read_text() == "chr22\t100\tA\tG\nchrX\t5000000\tC\tT\n"


def test_parse_mpileup_picks_the_site_alt():
    sites = {("chr22", 100): ("A", "G"), ("chr22", 200): ("C", "T")}
    text = "chr22\t100\tA\tT,G,<*>\t30,2,18,0\nchr22\t200\tC\t<*>\t40,0\nchr22\t999\tA\tG\t1,1\n"
    assert parse_mpileup(text, sites) == [("chr22", 100, 30, 18), ("chr22", 200, 40, 0)]


def test_infer_sex():
    assert infer_sex(het=0.01, y_ratio=0.4) == "M"
    assert infer_sex(het=0.3, y_ratio=0.01) == "F"
    assert infer_sex(het=0.3, y_ratio=0.4) == "?"


def _write_sample(path, genotypes, y, rng, x_het):
    with open(path, "w") as fo:
        fo.write(f"#depth\tchrY={y}\tchr22=100.0\n")
        fo.write("chrom\tpos\tref_n\talt_n\n")
        for pos, g in enumerate(genotypes):
            alt = sum(rng.random() < g / 2 for _ in range(60))
            fo.write(f"chr22\t{pos}\t{60 - alt}\t{alt}\n")
        for pos in range(200):
            alt = 30 if (x_het and pos % 2) else 0
            fo.write(f"chrX\t{pos}\t{60 - alt}\t{alt}\n")


def test_analyze_flags_swapped_pair_and_sex(tmp_path):
    rng = random.Random(1)
    people = {p: [rng.choice((0, 1, 2)) for _ in range(500)] for p in "ABC"}
    info = [  # Br3's DLPFC tube actually holds person A
        ("D1", "Br1_DLPFC", "M", "A"), ("D2", "Br1_Caudate", "M", "A"),
        ("D3", "Br2_DLPFC", "F", "B"), ("D4", "Br2_Caudate", "F", "B"),
        ("D5", "Br3_DLPFC", "F", "A"), ("D6", "Br3_Caudate", "F", "C"),
    ]
    with open(tmp_path / "info.tsv", "w") as fo:
        fo.write("sample_id\tbr_region\tsex\n")
        for sid, br, sex, person in info:
            fo.write(f"{sid}\t{br}\t{sex}\n")
            male = person == "A"
            _write_sample(tmp_path / f"{sid}.tsv", people[person], 40.0 if male else 0.5, rng, not male)

    assert main(["analyze", "--sample-info", str(tmp_path / "info.tsv"), "--out-dir", str(tmp_path)]) == 0
    pairs = {r["sample"]: r for r in csv.DictReader(open(tmp_path / "identity_pairs.tsv"), delimiter="\t")}
    assert pairs["D3"]["status"] == "OK" and pairs["D4"]["status"] == "OK"
    assert pairs["D5"]["status"] == "SWAP?" and pairs["D5"]["best_match"] in ("Br1_DLPFC", "Br1_Caudate")
    assert pairs["D6"]["status"] == "PAIR_MISMATCH"
    assert pairs["D1"]["status"] == "DUPLICATE_OR_RELATIVE?"   # D5 is also person A
    sex = {r["sample"]: r for r in csv.DictReader(open(tmp_path / "identity_sex.tsv"), delimiter="\t")}
    assert sex["D1"]["status"] == "OK" and sex["D3"]["status"] == "OK"
    assert sex["D5"]["status"] == "MISMATCH" and sex["D5"]["sex_inferred"] == "M"
