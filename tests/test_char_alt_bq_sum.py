"""Characterization test for scripts/alt_bq_sum.py against M1b golden.

Activates once ``tests/data/golden/filtering/alt_bq_sum/`` exists.
See test_char_somatic_vaf_filter.py for the contract.

DEF-002 of the v1.1.0 audit reclassified utils/alt_bq_sum.py from MOVE
to MERGE_INTO because it IS output-affecting via the BQ pileup path
when ``load_config`` is invoked. The strip is provably output-preserving
when ``--reference`` / ``--conda-env`` are not supplied (the documented
mode for tests/data/) — this golden test pins that contract.
"""

from __future__ import annotations

from pathlib import Path

import pytest

GOLDEN_DIR = Path(__file__).parent / "data" / "golden" / "filtering" / "alt_bq_sum"

pytestmark = pytest.mark.skipif(
    not GOLDEN_DIR.exists(),
    reason="M1b golden snapshots not yet generated — run on cluster",
)


def test_alt_bq_sum_matches_golden(tmp_path: Path) -> None:
    raise AssertionError(
        "M1b goldens present but characterization test body not yet authored"
    )
