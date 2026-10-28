"""Tests for scripts/indel_census.py (no bcftools needed)."""

from scripts.indel_census import parse_line, stage_of


def rec(**kw):
    base = "chr1\t100\tA\tAT\t1\tT\t3,4\t5.2\t0.08\t90,8"
    r = parse_line(base)
    r.update(kw)
    return r


def test_parse_line_fields():
    r = parse_line("chr1\t100\tACG\tA\t.\t.\t.\t.\t0.12\t70,9\n")
    assert (r["length"], r["alt_reads"], r["is_str"], r["popaf_v"]) == (-2, 9, False, 10.0)
    assert parse_line("chr1\t100\tA\tAT\t1\tT\t3,4\t5.2\t0.08\t90,8")["is_str"] is True


def test_parse_line_skips_multiallelic_and_malformed():
    assert parse_line("chr1\t100\tA\tAT,ATT\t.\t.\t.\t.\t0.1,0.1\t80,5,5") is None
    assert parse_line("chr1\t100\tA") is None


def test_stage_cascade():
    assert stage_of(rec(), accessible=False) == ["pass"]
    assert stage_of(rec(), accessible=True) == ["pass", "accessible", "rare", "low_vaf"]
    assert stage_of(rec(popaf_v=2.0), accessible=True) == ["pass", "accessible"]       # gnomAD AF 1%
    assert stage_of(rec(vaf=0.45), accessible=True) == ["pass", "accessible", "rare"]   # germline het
    assert stage_of(rec(alt_reads=3), accessible=True) == ["pass", "accessible", "rare"]
