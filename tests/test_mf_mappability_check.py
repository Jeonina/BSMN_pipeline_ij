"""Tests for scripts/mf_mappability_check.py pure logic (no apptainer)."""

from scripts.mf_mappability_check import (
    bed_line,
    parse_bw_average,
    parse_predictions,
    read_final,
    swap_mappability,
)


def test_bed_line_matches_mf_rule():
    # ReadLevel_Features_extraction.py: chrom, pos, pos+1, id; 'chr' added when missing
    assert bed_line("D1~chr1~100~A~G") == "chr1\t100\t101\tD1~chr1~100~A~G"
    assert bed_line("D1~7~5~C~T") == "chr7\t5\t6\tD1~7~5~C~T"


def test_parse_bw_average_reads_column_6():
    text = "a\t1\t0\t0\t0\t0\nb\t1\t1\t0.5\t0.5\t0.958333\n"
    assert parse_bw_average(text) == {"a": 0.0, "b": 0.958333}


def test_parse_predictions_label_and_prob_skip_header():
    row = ["D1~chr1~100~A~G"] + ["x"] * 33 + ["mosaic", "x", "0.81"]
    text = "\t".join(["id"] + ["h"] * 36) + "\n" + "\t".join(f'"{v}"' if i == 0 else v for i, v in enumerate(row)) + "\n"
    assert parse_predictions(text) == {"D1~chr1~100~A~G": ("mosaic", "0.81")}


def test_swap_mappability_replaces_column_3_only(tmp_path):
    lines = ["id\tconflict_num\tmappability\ttype\n", "D1~chr1~100~A~G\t0\t0.0\tSNP\n", "D1~chr2~5~C~T\t0\t1.0\tSNP\n"]
    out = swap_mappability(lines, {"D1~chr1~100~A~G": 0.958333})
    assert out == [lines[0], "D1~chr1~100~A~G\t0\t0.958333\tSNP\n", lines[2]]


def test_read_final_keys_by_sample(tmp_path):
    f = tmp_path / "D1.final.txt"
    f.write_text("chr1\t100\tA\tG\n\n")                       # final.txt has exactly 4 columns
    assert read_final(str(f), "D1") == {"D1~chr1~100~A~G"}
    assert read_final(str(tmp_path / "missing.txt"), "D1") == set()
