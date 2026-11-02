"""Tests for scripts/phase_somatic.py input building (no container needed)."""

from scripts.phase_somatic import position_lines, read_final


def test_read_final_keeps_snvs_only(tmp_path):
    f = tmp_path / "S.final.txt"
    f.write_text("chr1\t100\tA\tG\textra\nchr1\t200\tAT\tA\nchr2\t50\tC\tT\n#comment\n")
    assert read_final(str(f)) == [("chr1", 100, "A", "G"), ("chr2", 50, "C", "T")]


def test_position_lines_bed_like_and_sorted():
    calls = {"S2": [("chr2", 50, "C", "T")], "S1": [("chr1", 300, "G", "A"), ("chr1", 100, "A", "G")]}
    assert position_lines(calls) == [
        "chr1\t99\t100\tA\tG\tS1",
        "chr1\t299\t300\tG\tA\tS1",
        "chr2\t49\t50\tC\tT\tS2",
    ]
