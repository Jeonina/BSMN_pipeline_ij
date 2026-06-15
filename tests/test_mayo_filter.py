"""Unit tests for scripts/mayo_filter.py decision logic.

Covers the four BSMN E.mayo_filters.sh criteria (repeat, multiallelic,
both-strand, strand-balance) plus the TSV join. No BAM/CRAM needed — the
strand/repeat computations live in the upstream scripts.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).parent.parent / "scripts"
_spec = importlib.util.spec_from_file_location("mayo_filter", SCRIPTS_DIR / "mayo_filter.py")
assert _spec and _spec.loader
mayo = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mayo)


def _strand(**over):
    """A clean, passing strand record; override fields per test."""
    rec = {
        "chr": "chr20", "pos": "47053475", "ref": "G", "alt": "T",
        "total": 273, "total_fwd": 147, "total_rev": 126,
        "ref_n": 242, "ref_fwd": 130, "ref_rev": 112,
        "alt_n": 31, "alt_fwd": 17, "alt_rev": 14,
        "p_poisson": 0.30, "p_fisher": 0.56,
    }
    rec.update(over)
    return rec


def _repeat(**over):
    rec = {"chr": "chr20", "pos": "47053475", "ref": "G", "alt": "T",
           "repeat_n": 1, "repeat_length": 1, "repeat_seq": "[G>T]"}
    rec.update(over)
    return rec


# --- passes_mayo: one criterion at a time ----------------------------------

def test_clean_record_passes():
    kept, reason = mayo.passes_mayo(_strand(), _repeat(), 4, 10, 0.05)
    assert kept is True
    assert reason is None


def test_repeat_n_at_threshold_rejected():
    kept, reason = mayo.passes_mayo(_strand(), _repeat(repeat_n=4), 4, 10, 0.05)
    assert (kept, reason) == (False, "repeat")


def test_repeat_length_at_threshold_rejected():
    kept, reason = mayo.passes_mayo(_strand(), _repeat(repeat_length=10), 4, 10, 0.05)
    assert (kept, reason) == (False, "repeat")


def test_multiallelic_rejected():
    # total != ref_n + alt_n  -> a third allele is present
    kept, reason = mayo.passes_mayo(_strand(total=300), _repeat(), 4, 10, 0.05)
    assert (kept, reason) == (False, "multiallelic")


def test_single_strand_alt_rejected():
    kept, reason = mayo.passes_mayo(
        _strand(alt_fwd=31, alt_rev=0, ref_n=242, total=273), _repeat(), 4, 10, 0.05
    )
    assert (kept, reason) == (False, "single_strand")


def test_strand_bias_rejected_when_both_p_below():
    kept, reason = mayo.passes_mayo(
        _strand(p_poisson=0.01, p_fisher=0.02), _repeat(), 4, 10, 0.05
    )
    assert (kept, reason) == (False, "strand_bias")


def test_strand_balance_passes_if_either_p_high():
    # poisson low but fisher high -> kept (OR semantics)
    kept, _ = mayo.passes_mayo(_strand(p_poisson=0.001, p_fisher=0.10), _repeat(), 4, 10, 0.05)
    assert kept is True


def test_criterion_order_repeat_reported_first():
    # repeat AND strand_bias both fail -> repeat reported (evaluated first)
    kept, reason = mayo.passes_mayo(
        _strand(p_poisson=0.0, p_fisher=0.0), _repeat(repeat_n=9), 4, 10, 0.05
    )
    assert (kept, reason) == (False, "repeat")


# --- read_tsv: header parsing + typing --------------------------------------

def test_read_tsv_join_and_types(tmp_path):
    strand = tmp_path / "s.tsv"
    strand.write_text(
        "#chr\tpos\tref\talt\ttotal\ttotal_fwd\ttotal_rev\ttotal_ratio\tp_poisson\t"
        "ref_n\tref_fwd\tref_rev\tref_ratio\talt_n\talt_fwd\talt_rev\talt_ratio\tp_fisher\n"
        "chr20\t47053475\tG\tT\t273\t147\t126\t1.17\t0.30\t"
        "242\t130\t112\t1.16\t31\t17\t14\t1.21\t0.56\n"
    )
    repeat = tmp_path / "r.tsv"
    repeat.write_text(
        "#chr\tpos\tref\talt\trepeat_n\trepeat_length\trepeat_seq\n"
        "chr20\t47053475\tG\tT\t1\t1\t[G>T]\n"
    )
    s = mayo.read_tsv(str(strand))
    r = mayo.read_tsv(str(repeat))
    key = ("chr20", "47053475", "G", "T")
    assert s[key]["total"] == 273 and isinstance(s[key]["total"], int)
    assert s[key]["p_fisher"] == pytest.approx(0.56)
    assert r[key]["repeat_n"] == 1
    kept, reason = mayo.passes_mayo(s[key], r[key], 4, 10, 0.05)
    assert kept is True and reason is None
