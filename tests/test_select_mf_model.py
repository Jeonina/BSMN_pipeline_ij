"""Unit tests for scripts/select_mf_model.py pure logic.

Covers model discovery, nearest-depth choice, samtools coverage parsing and the
auto/pinned decision. The `samtools coverage` container call is not tested.
"""

from __future__ import annotations

import importlib.util
import logging
from pathlib import Path

import pytest

SCRIPTS_DIR = Path(__file__).parent.parent / "scripts"
_spec = importlib.util.spec_from_file_location(
    "select_mf_model", SCRIPTS_DIR / "select_mf_model.py"
)
assert _spec and _spec.loader
sm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sm)


@pytest.fixture
def models_dir(tmp_path):
    # The file set parklab/MosaicForecast ships in models_trained/.
    for name in [
        "50xRFmodel_addRMSK_Refine.rds",
        "100xRFmodel_addRMSK_Refine.rds",
        "150xRFmodel_addRMSK_Refine.rds",
        "200xRFmodel_addRMSK_Refine.rds",
        "250xRFmodel_addRMSK_Refine.rds",
        "brain_MT2-PON.250x_MosaicForecast-Phase_pvalue.rds",
        "README.md",
    ]:
        (tmp_path / name).touch()
    return tmp_path


# --- list_models -------------------------------------------------------------


def test_list_models_refine_only(models_dir):
    models = sm.list_models(str(models_dir), "Refine")
    assert sorted(models) == [50, 100, 150, 200, 250]
    assert models[100].endswith("100xRFmodel_addRMSK_Refine.rds")


def test_list_models_excludes_other_families(models_dir):
    # The brain Phase-pvalue model does not pair with mode Refine.
    models = sm.list_models(str(models_dir), "Refine")
    assert not any("brain" in p for p in models.values())


def test_list_models_no_match_for_other_mode(models_dir):
    assert sm.list_models(str(models_dir), "Phase") == {}


# --- nearest_depth -----------------------------------------------------------


@pytest.mark.parametrize(
    ("depth", "expected"),
    [(98.99, 100), (110.5, 100), (30.0, 50), (172.0, 150), (400.0, 250), (190.0, 200)],
)
def test_nearest_depth(depth, expected):
    assert sm.nearest_depth(depth, [50, 100, 150, 200, 250]) == expected


def test_nearest_depth_tie_goes_lower():
    assert sm.nearest_depth(125.0, [150, 100]) == 100


def test_nearest_depth_empty():
    with pytest.raises(ValueError):
        sm.nearest_depth(100.0, [])


# --- parse_meandepth ---------------------------------------------------------


def test_parse_meandepth():
    text = (
        "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
        "chr20\t1\t64444167\t42000000\t63000000\t97.8\t98.99\t36.1\t58.2\n"
    )
    assert sm.parse_meandepth(text) == pytest.approx(98.99)


def test_parse_meandepth_no_data_row():
    with pytest.raises(ValueError):
        sm.parse_meandepth("#rname\tmeandepth\n")


# --- choose ------------------------------------------------------------------


def test_choose_auto(models_dir):
    models = sm.list_models(str(models_dir), "Refine")
    model, label, source = sm.choose(58.0, models, "auto")
    assert (label, source) == (50, "auto")
    assert model.endswith("50xRFmodel_addRMSK_Refine.rds")


def test_choose_auto_without_models():
    with pytest.raises(ValueError):
        sm.choose(100.0, {}, "auto")


def test_choose_pinned_matching_depth_no_warning(models_dir, caplog):
    models = sm.list_models(str(models_dir), "Refine")
    pinned = str(models_dir / "100xRFmodel_addRMSK_Refine.rds")
    with caplog.at_level(logging.WARNING):
        assert sm.choose(104.0, models, pinned) == (pinned, 100, "pinned")
    assert not caplog.records


def test_choose_pinned_mismatch_warns_but_keeps_pin(models_dir, caplog):
    models = sm.list_models(str(models_dir), "Refine")
    pinned = str(models_dir / "250xRFmodel_addRMSK_Refine.rds")
    with caplog.at_level(logging.WARNING):
        assert sm.choose(99.0, models, pinned) == (pinned, 250, "pinned")
    assert "nearest 100x" in caplog.text
