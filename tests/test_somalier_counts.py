"""Tests for scripts/somalier_counts.py (no bcftools needed)."""

import gzip
import struct

from scripts.somalier_counts import load_sites, parse_counts, write_somalier

SITES = [("chr1", 100, "G", "C"), ("chr2", 200, "T", "C"), ("chrX", 300, "A", "G"), ("chrY", 400, "C", "T")]


def test_load_sites_keeps_vcf_order(tmp_path):
    vcf = tmp_path / "s.vcf.gz"
    with gzip.open(vcf, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for c, p, r, a in SITES:
            fh.write(f"{c}\t{p}\t.\t{r}\t{a}\t.\tPASS\tAF=0.3\n")
    assert load_sites(str(vcf)) == SITES


def test_parse_counts_splits_ref_alt_other():
    text = "chr1\t100\tG\tA,C,<*>\t30,2,18,0\nchr2\t200\tT\t<*>\t40,0\nchr9\t1\tA\tG\t1,1\n"
    assert parse_counts(text, SITES) == {("chr1", 100): (30, 18, 2), ("chr2", 200): (40, 0, 0)}


def test_write_somalier_binary_layout(tmp_path):
    out = tmp_path / "D1.somalier"
    write_somalier(str(out), "D1", SITES, {("chr1", 100): (30, 18, 2), ("chrX", 300): (5, 6, 0)})
    b = out.read_bytes()
    assert b[:4] == bytes([2, 2]) + b"D1"
    assert struct.unpack("<HHH", b[4:10]) == (2, 1, 1)          # auto, X, Y
    recs = [struct.unpack("<III", b[10 + 12 * i: 22 + 12 * i]) for i in range(4)]
    assert recs == [(30, 18, 2), (0, 0, 0), (5, 6, 0), (0, 0, 0)]
    assert len(b) == 2 + 2 + 6 + 12 * 4
