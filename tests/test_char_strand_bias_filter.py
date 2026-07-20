"""Characterization test for scripts/strand_bias_filter.py against M1b golden.

Activates once ``tests/data/golden/filtering/strand_bias/`` exists.
See test_char_somatic_vaf_filter.py for the contract.
"""

from __future__ import annotations

from pathlib import Path

import pytest

GOLDEN_DIR = Path(__file__).parent / "data" / "golden" / "filtering" / "strand_bias"

pytestmark = pytest.mark.skipif(
    not GOLDEN_DIR.exists(),
    reason="M1b golden snapshots not yet generated — run on cluster",
)


def test_strand_bias_filter_matches_golden(tmp_path: Path) -> None:
    raise AssertionError("M1b goldens present but characterization test body not yet authored")
