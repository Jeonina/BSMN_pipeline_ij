"""Unit tests for scripts/cnvnator_filter.py pure logic.

Covers region construction and genotype-output parsing (the bug-prone parts).
The `cnvnator -genotype` container call needs a ROOT file and is not tested.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

SCRIPTS_DIR = Path(__file__).parent.parent / "scripts"
_spec = importlib.util.spec_from_file_location(
    "cnvnator_filter", SCRIPTS_DIR / "cnvnator_filter.py"
)
assert _spec and _spec.loader
cf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cf)


# --- region_str -------------------------------------------------------------


def test_region_window():
    assert cf.region_str("chr20", 47053475, 1000) == "chr20:47052475-47054475"


def test_region_start_clamped_to_one():
    assert cf.region_str("chr1", 500, 1000) == "chr1:1-1500"


# --- parse_genotype ---------------------------------------------------------


def test_parse_genotype_default_field():
    # "Genotype <region> <normRD> <CN> ..." → CN at index 3 (BSMN $9)
    text = (
        "Genotype chr20:47052475-47054475 1.01 2.0 0.5\n"
        "Genotype chr20:5121314-5123314 1.48 3.1 0.2\n"
    )
    got = cf.parse_genotype(text)
    assert got["chr20:47052475-47054475"] == 2.0
    assert got["chr20:5121314-5123314"] == 3.1


def test_parse_genotype_ignores_noise_and_bad_lines():
    text = (
        "Assuming male\n"
        "Genotype chr20:1-2000 0.9 2.0\n"
        "Genotype chr20:9-9 onlytwo\n"  # too short → skipped
        "random log line\n"
    )
    got = cf.parse_genotype(text)
    assert got == {"chr20:1-2000": 2.0}


def test_parse_genotype_custom_field():
    text = "Genotype chrX:100-200 5.5 9.9\n"
    assert cf.parse_genotype(text, cn_field=2) == {"chrX:100-200": 5.5}


def test_read_candidates(tmp_path):
    p = tmp_path / "c.txt"
    p.write_text("chr20\t47053475\tG\tT\nchr20\t5122314\ta\tt\n")
    assert cf.read_candidates(str(p)) == [
        ("chr20", 47053475, "G", "T"),
        ("chr20", 5122314, "A", "T"),
    ]
