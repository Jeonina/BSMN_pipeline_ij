"""Tests for scripts/cohort_audit.py (cohort-wide once-only completeness audit)."""

import csv
import gzip
import os

from scripts.cohort_audit import main

FLAGSTAT = """\
1000 + 0 in total (QC-passed reads + QC-failed reads)
960 + 0 primary
96 + 0 primary duplicates
940 + 0 primary mapped (97.92% : N/A)
"""


def _w(path, text=""):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def _sample(root, s, host_logs, final_lines=2):
    m = root / "results" / "mapping" / s
    _w(str(m / f"{s}.cram"))
    _w(str(m / f"{s}.cram.crai"))
    _w(str(m / "flagstat.txt"), FLAGSTAT)
    c = root / "results" / "calling" / s
    os.makedirs(c)
    with gzip.open(c / f"{s}.filtered.vcf.gz", "wt") as fh:
        fh.write(
            f"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{s}\n"
        )
    _w(str(c / f"{s}.filtered.vcf.gz.tbi"))
    _w(
        str(root / "results" / "filtering" / s / f"{s}.final.txt"),
        "".join(f"chr1\t{i}\tA\tC\n" for i in range(final_lines)),
    )
    for logs in host_logs:
        _w(str(logs / "calling" / s / "filter_mutect_calls.log"))
        _w(
            str(logs / "calling" / s / "calculate_contamination.log"),
            "contamination_value=0.0012\n",
        )
        f = logs / "filtering" / s
        _w(str(f / "accessibility_filter.log"), "variants_input=500\nvariants_kept=400\n")
        _w(str(f / "germline_filter.log"), "variants_kept=100\n")
        _w(str(f / "vaf_filter.log"), "variants_kept=20\n")
        _w(str(f / "cnvnator_filter.log"), "cnvnator: 20 -> 15  (removed)\n")
        _w(
            str(f / "mosaicforecast_filter.log"),
            f"candidates in=15  sample={s}  model=/m/250x.rds\n"
            "mosaicforecast: 15 -> 3 (min_prob=0.5)\n",
        )
        _w(str(f / "pon_mask_filter.log"), "variants_kept=2\nEND pon_mask_filter\n")


def _run(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    out = tmp_path / "audit.tsv"
    rc = main(
        [
            "--barcode-map",
            "map.tsv",
            "--results",
            "results",
            "--log-dir",
            "vm100=l100",
            "--log-dir",
            "vm102=l102",
            "--out",
            str(out),
        ]
    )
    with open(out) as fh:
        return rc, {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}


def test_clean_cohort_split_across_hosts(tmp_path, monkeypatch):
    _w(str(tmp_path / "map.tsv"), "sample_id\tfc\tudb\nD1\tF1\tU1\nD1\tF1\tU2\nD2\tF2\tU1\n")
    _sample(tmp_path, "D1", [tmp_path / "l100"])
    _sample(tmp_path, "D2", [tmp_path / "l102"])
    rc, rows = _run(tmp_path, monkeypatch)
    assert rc == 0
    assert rows["D1"]["filtering_host"] == "vm100" and rows["D2"]["filtering_host"] == "vm102"
    d1 = rows["D1"]
    assert (d1["pass_snv"], d1["vaf"], d1["mosaicforecast"], d1["pon"], d1["final_snv"]) == (
        "500",
        "20",
        "3",
        "2",
        "2",
    )
    assert d1["mf_model"] == "250x.rds" and d1["contamination"] == "0.0012"
    assert d1["n_rg_expected"] == "2" and d1["dup_pct"] == "10.0"


def test_cascade_falls_back_to_other_host_logs(tmp_path, monkeypatch):
    _w(str(tmp_path / "map.tsv"), "sample_id\tfc\tudb\nD1\tF1\tU1\n")
    _sample(tmp_path, "D1", [tmp_path / "l102"])  # pilot-style: early logs elsewhere
    _w(
        str(tmp_path / "l100" / "filtering" / "D1" / "pon_mask_filter.log"),
        "variants_kept=2\nEND pon_mask_filter\n",
    )
    os.remove(tmp_path / "l102" / "filtering" / "D1" / "pon_mask_filter.log")
    rc, rows = _run(tmp_path, monkeypatch)
    assert rows["D1"]["filtering_host"] == "vm100"
    assert rows["D1"]["pass_snv"] == "500" and rows["D1"]["pon"] == "2"


def test_flags_duplicate_missing_and_mismatch(tmp_path, monkeypatch):
    _w(str(tmp_path / "map.tsv"), "sample_id\tfc\tudb\nD1\tF1\tU1\nD2\tF2\tU1\nD3\tF3\tU1\n")
    _sample(tmp_path, "D1", [tmp_path / "l100", tmp_path / "l102"])  # ran twice
    _sample(tmp_path, "D2", [tmp_path / "l100"], final_lines=5)  # final != pon log
    rc, rows = _run(tmp_path, monkeypatch)
    assert rc == 1
    assert "filtering_ran_on_vm100+vm102" in rows["D1"]["problems"]
    assert "final_lines(5)!=pon_log(2)" in rows["D2"]["problems"]
    assert "mapping_incomplete" in rows["D3"]["problems"]
