"""Unit tests for scripts/merge_donor_crams.py (pure helpers; no samtools needed)."""

import pytest

from scripts.merge_donor_crams import load_donors, rewrite_sm


def _info(tmp_path, rows):
    p = tmp_path / "sample_info.tsv"
    p.write_text("sample_id\tbr_region\n" + "".join(f"{s}\t{b}\n" for s, b in rows))
    return str(p)


def test_load_donors_pairs_regions(tmp_path):
    path = _info(tmp_path, [("D1", "Br0846_DLPFC"), ("D2", "Br0846_Caudate")])
    assert load_donors(path) == {"Br0846": {"DLPFC": "D1", "Caudate": "D2"}}


def test_load_donors_rejects_missing_region(tmp_path):
    path = _info(tmp_path, [("D1", "Br0846_DLPFC")])
    with pytest.raises(ValueError, match="without both regions"):
        load_donors(path)


def test_load_donors_rejects_duplicate_region(tmp_path):
    path = _info(tmp_path, [("D1", "Br0846_DLPFC"), ("D2", "Br0846_DLPFC")])
    with pytest.raises(ValueError, match="two DLPFC"):
        load_donors(path)


def test_load_donors_rejects_unknown_region(tmp_path):
    path = _info(tmp_path, [("D1", "Br0846_Hippo")])
    with pytest.raises(ValueError, match="unexpected br_region"):
        load_donors(path)


HEADER = (
    "@HD\tVN:1.6\tSO:coordinate\n"
    "@SQ\tSN:chr1\tLN:248956422\n"
    "@RG\tID:D1.FC1_UDB-1\tSM:D1\tPL:illumina\tLB:D1\tPU:FC1_UDB-1\n"
    "@RG\tID:D2.FC2_UDB-2\tSM:D2\tPL:illumina\tLB:D2\tPU:FC2_UDB-2\n"
    "@PG\tID:samtools\tPN:samtools\n"
)


def test_rewrite_sm_sets_donor_and_keeps_rest():
    out = rewrite_sm(HEADER, "Br0846", {"D1", "D2"})
    rg = [line for line in out.splitlines() if line.startswith("@RG")]
    assert all("\tSM:Br0846\t" in line for line in rg)
    assert "ID:D1.FC1_UDB-1" in rg[0] and "LB:D1" in rg[0]   # ID/LB untouched
    assert out.replace("SM:Br0846", "SM:X") == HEADER.replace("SM:D1", "SM:X").replace("SM:D2", "SM:X")


def test_rewrite_sm_rejects_foreign_sample():
    with pytest.raises(ValueError, match="SM set"):
        rewrite_sm(HEADER, "Br0846", {"D1", "D3"})
