"""Tests for scripts.validate_resolved_params (REQ-SNK-009 / AC-18).

The resolved_params.yaml schema is documented in PIPELINE.md §169-177 and
emitted by scripts/auto_params.py. The six required top-level keys are:

  - resolved_at                       (str, ISO-8601)
  - input_fastq                       (str, absolute path)
  - sequencer                         (str)
  - bwa_threads                       (int, >= 1)
  - optical_duplicate_pixel_distance  (int, >= 0)
  - bqsr_memory_gb                    (int, >= 1)

Additional keys (e.g. ``sequencer_evidence``, ``sort_threads``) MAY be
present; the validator MUST NOT reject them.
"""

from __future__ import annotations

from pathlib import Path

import pytest

# pyyaml is required at validator-import time; mark module as needing it
yaml = pytest.importorskip("yaml")

from scripts import validate_resolved_params as v  # noqa: E402

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _minimal_valid() -> dict[str, object]:
    return {
        "resolved_at": "2026-05-12T10:00:00",
        "input_fastq": "/abs/path/sample_R1.fastq.gz",
        "sequencer": "NovaSeq 6000",
        "bwa_threads": 8,
        "optical_duplicate_pixel_distance": 2500,
        "bqsr_memory_gb": 16,
    }


def _write_yaml(tmp_path: Path, payload: object) -> Path:
    p = tmp_path / "resolved_params.yaml"
    p.write_text(yaml.safe_dump(payload), encoding="utf-8")
    return p


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------


def test_valid_minimal_yaml(tmp_path: Path) -> None:
    p = _write_yaml(tmp_path, _minimal_valid())
    out = v.validate(p)
    assert out.sequencer == "NovaSeq 6000"
    assert out.bwa_threads == 8
    assert out.optical_duplicate_pixel_distance == 2500


def test_extra_keys_are_allowed(tmp_path: Path) -> None:
    payload = _minimal_valid()
    payload["sequencer_evidence"] = {"source": "instrument_id"}
    payload["sort_threads"] = 4
    p = _write_yaml(tmp_path, payload)
    out = v.validate(p)
    assert out.sequencer == "NovaSeq 6000"


def test_optical_zero_is_allowed(tmp_path: Path) -> None:
    payload = _minimal_valid()
    payload["optical_duplicate_pixel_distance"] = 0
    p = _write_yaml(tmp_path, payload)
    out = v.validate(p)
    assert out.optical_duplicate_pixel_distance == 0


# ---------------------------------------------------------------------------
# Failure paths — each must raise SystemExit(2)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "missing_key",
    [
        "resolved_at",
        "input_fastq",
        "sequencer",
        "bwa_threads",
        "optical_duplicate_pixel_distance",
        "bqsr_memory_gb",
    ],
)
def test_missing_required_key_exits_2(tmp_path: Path, missing_key: str) -> None:
    payload = _minimal_valid()
    payload.pop(missing_key)
    p = _write_yaml(tmp_path, payload)
    with pytest.raises(SystemExit) as exc:
        v.validate(p)
    assert exc.value.code == 2


def test_wrong_type_str_field_exits_2(tmp_path: Path) -> None:
    payload = _minimal_valid()
    payload["resolved_at"] = 12345  # int instead of str
    p = _write_yaml(tmp_path, payload)
    with pytest.raises(SystemExit) as exc:
        v.validate(p)
    assert exc.value.code == 2


def test_wrong_type_int_field_exits_2(tmp_path: Path) -> None:
    payload = _minimal_valid()
    payload["bwa_threads"] = "eight"
    p = _write_yaml(tmp_path, payload)
    with pytest.raises(SystemExit) as exc:
        v.validate(p)
    assert exc.value.code == 2


def test_negative_bwa_threads_exits_2(tmp_path: Path) -> None:
    payload = _minimal_valid()
    payload["bwa_threads"] = 0
    p = _write_yaml(tmp_path, payload)
    with pytest.raises(SystemExit) as exc:
        v.validate(p)
    assert exc.value.code == 2


def test_negative_optical_exits_2(tmp_path: Path) -> None:
    payload = _minimal_valid()
    payload["optical_duplicate_pixel_distance"] = -1
    p = _write_yaml(tmp_path, payload)
    with pytest.raises(SystemExit) as exc:
        v.validate(p)
    assert exc.value.code == 2


def test_malformed_yaml_exits_2(tmp_path: Path) -> None:
    p = tmp_path / "bad.yaml"
    p.write_text("resolved_at: : :\n  not - valid", encoding="utf-8")
    with pytest.raises(SystemExit) as exc:
        v.validate(p)
    assert exc.value.code == 2


def test_missing_file_exits_2(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as exc:
        v.validate(tmp_path / "does_not_exist.yaml")
    assert exc.value.code == 2


def test_top_level_not_mapping_exits_2(tmp_path: Path) -> None:
    p = _write_yaml(tmp_path, ["not", "a", "mapping"])
    with pytest.raises(SystemExit) as exc:
        v.validate(p)
    assert exc.value.code == 2


def test_bool_is_rejected_as_int(tmp_path: Path) -> None:
    # In Python, bool is a subclass of int — guard against True/False slipping
    # into bwa_threads.
    payload = _minimal_valid()
    payload["bwa_threads"] = True  # type: ignore[assignment]
    p = _write_yaml(tmp_path, payload)
    with pytest.raises(SystemExit) as exc:
        v.validate(p)
    assert exc.value.code == 2
