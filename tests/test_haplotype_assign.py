"""Tests for scripts/haplotype_assign.py (no pysam / CRAM needed)."""

from collections import Counter

from scripts.haplotype_assign import classify, count_pairs, donor_column, pair_allele, pair_haplotype


def test_pair_haplotype_needs_agreeing_votes():
    assert pair_haplotype(["A", "A"]) == "A"
    assert pair_haplotype(["A", "B"]) is None
    assert pair_haplotype([]) is None


def test_pair_allele_drops_other_bases_and_disagreeing_mates():
    assert pair_allele({"T"}, "C", "T") == "alt"
    assert pair_allele({"C"}, "C", "T") == "ref"
    assert pair_allele({"C", "T"}, "C", "T") is None
    assert pair_allele({"G"}, "C", "T") is None


def test_classify_follows_bae2022_definitions():
    def n(**kw):
        return Counter({tuple(k.split("_")): v for k, v in kw.items()})
    assert classify(n(alt_A=3, ref_A=20, ref_B=25), 2) == "supported"
    assert classify(n(alt_A=3, alt_B=4, ref_A=10, ref_B=9), 2) == "unsupported"
    assert classify(n(alt_A=3, alt_B=1, ref_A=10), 2) == "ambiguous"
    assert classify(n(alt_B=12, ref_A=11), 2) == "germline_like"
    assert classify(n(ref_A=10, ref_B=10), 2) == "unphased"


def test_donor_column_matches_donor_prefix():
    assert donor_column(["Br5094_Caudate", "Br5458_DLPFC"], "Br5458") == "Br5458_DLPFC"
    assert donor_column(["Br5094_Caudate"], "Br509") is None


class _Read:
    def __init__(self, name, start, seq, mapq=60):
        self.query_name, self._start, self.query_sequence = name, start, seq
        self.query_qualities = [30] * len(seq)
        self.mapping_quality = mapq
        self.is_unmapped = self.is_secondary = self.is_supplementary = self.is_duplicate = self.is_qcfail = False

    def get_aligned_pairs(self, matches_only=True):
        return [(i, self._start + i) for i in range(len(self.query_sequence))]


class _Bam:
    def __init__(self, reads):
        self.reads = reads

    def fetch(self, chrom, lo, hi):
        return iter(self.reads)


def test_count_pairs_uses_mates_and_phased_snps():
    # SNV at pos 11 (0-based 10) C>T; phased het SNP at pos 16: hap A 'G', hap B 'A'
    snps = {16: ("G", "A")}
    reads = [_Read("p1", 0, "AAAAAAAAAATAAAAGAAAA"),   # alt, hap A
             _Read("p2", 0, "AAAAAAAAAACAAAAAAAAA"),   # ref, hap B
             _Read("p3", 0, "AAAAAAAAAATAAAAAAAAA"),   # alt, hap B
             _Read("p4", 0, "AAAAAAAAAATAAAA"),        # alt, no SNP covered -> unphased pair, dropped
             _Read("p5", 0, "AAAAAAAAAACAAAAGAAAA"), _Read("p5", 12, "AAAAAAAA")]   # mates disagree at SNP -> dropped
    n = count_pairs(_Bam(reads), "chr1", 11, "C", "T", snps, 20, 20)
    assert n == Counter({("alt", "A"): 1, ("ref", "B"): 1, ("alt", "B"): 1})
