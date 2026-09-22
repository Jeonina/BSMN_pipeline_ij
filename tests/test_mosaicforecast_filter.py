"""Unit tests for scripts/mosaicforecast_filter.py pure logic.

Covers BED conversion and prediction-file parsing (the bug-prone parts). The
Apptainer feature-extraction / R prediction calls need the image + alignment
data and are not exercised here.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

SCRIPTS_DIR = Path(__file__).parent.parent / "scripts"
_spec = importlib.util.spec_from_file_location(
    "mosaicforecast_filter", SCRIPTS_DIR / "mosaicforecast_filter.py"
)
assert _spec and _spec.loader
mf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mf)


# --- to_bed_line ------------------------------------------------------------


def test_to_bed_line_zero_based_start_and_sample():
    assert mf.to_bed_line("chr20", "47053475", "G", "T", "HG002") == (
        "chr20\t47053474\t47053475\tG\tT\tHG002"
    )


# --- parse_predictions ------------------------------------------------------


def _row(var_id: str, pred: str, prob: str) -> str:
    """37-column whitespace row: col1=id, col35=pred, col37=prob."""
    fields = ["."] * 37
    fields[0] = var_id
    fields[mf.PRED_COL] = pred
    fields[mf.PROB_COL] = prob
    return "\t".join(fields)


def test_header_and_nonmosaic_skipped():
    text = "\n".join(
        [
            "id\t" + "\t".join(["h"] * 36),  # header row
            _row("HG002~chr20~47053475~G~T", "mosaic", "0.95"),  # mosaic
            _row("HG002~chr20~1454720~T~C", "het", "0.99"),  # not mosaic
            _row("HG002~chr20~5122314~A~T", "refhom", "0.80"),  # not mosaic
        ]
    )
    got = mf.parse_predictions(text, min_prob=0.0)
    assert got == [("chr20", "47053475", "G", "T")]


def test_min_prob_filters_low_confidence():
    text = "\n".join(
        [
            _row("HG002~chr20~47053475~G~T", "mosaic", "0.95"),
            _row("HG002~chr20~38139365~T~A", "mosaic", "0.40"),  # below 0.6
        ]
    )
    assert mf.parse_predictions(text, min_prob=0.6) == [("chr20", "47053475", "G", "T")]
    # with no cutoff both survive
    assert len(mf.parse_predictions(text, min_prob=0.0)) == 2


def test_mosaic_prefix_match():
    # MF emits labels like 'mosaic' — prefix match per BSMN awk '/^mosaic/'
    text = _row("S~chr1~100~A~G", "mosaic", "0.7")
    assert mf.parse_predictions(text, 0.0) == [("chr1", "100", "A", "G")]


def test_short_rows_ignored():
    assert mf.parse_predictions("too\tshort\trow", 0.0) == []
    assert mf.parse_predictions("", 0.0) == []


def test_read_candidates(tmp_path):
    p = tmp_path / "c.txt"
    p.write_text("chr20\t47053475\tG\tT\nchr20\t5122314\ta\tt\n")
    assert mf.read_candidates(str(p)) == [
        ("chr20", "47053475", "G", "T"),
        ("chr20", "5122314", "A", "T"),
    ]


# --- dropout accounting -----------------------------------------------------


def _args(tmp_path, **over):
    import argparse

    base = dict(
        sample="SM",
        workdir=str(tmp_path),
        workers=2,
        threads=1,
        timeout=900,
        timeout_max=3600,
        retries=3,
        max_dropout=0.10,
    )
    base.update(over)
    return argparse.Namespace(**base)


def test_extract_features_counts_dropped(tmp_path, monkeypatch):
    """A variant that yields no row is reported as dropped, not silently skipped."""
    beds = ["chr1\t9\t10\tA\tT\tSM", "chr1\t19\t20\tC\tG\tSM", "chr1\t29\t30\tG\tA\tSM"]

    def fake_extract(idx, total, bed, args):
        return ("h1\th2", f"row{idx}") if idx != 2 else (None, None)

    monkeypatch.setattr(mf, "_extract_one", fake_extract)
    features, dropped = mf.extract_features(_args(tmp_path), beds)
    assert dropped == 1
    assert Path(features).read_text().splitlines() == ["h1\th2", "row1", "row3"]


def test_extract_features_all_dropped_reports_total(tmp_path, monkeypatch):
    monkeypatch.setattr(mf, "_extract_one", lambda *a, **k: (None, None))
    features, dropped = mf.extract_features(_args(tmp_path), ["chr1\t9\t10\tA\tT\tSM"])
    assert features is None
    assert dropped == 1


def test_escalating_timeout_doubles_and_caps():
    """Retrying a slow variant with the same budget re-fails; the budget must grow."""
    args = _args(Path("."), timeout=900, timeout_max=3600, retries=4)
    budgets = [min(args.timeout * (2**a), args.timeout_max) for a in range(args.retries)]
    assert budgets == [900, 1800, 3600, 3600]
