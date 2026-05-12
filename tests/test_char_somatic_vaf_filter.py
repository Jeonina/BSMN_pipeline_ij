"""Characterization test for scripts/somatic_vaf_filter.py against M1b golden.

Per SPEC-BSMN-REFACTOR-001 REQ-DUP-002 + REQ-TST-002, the M2 consolidation
strips ``load_config`` from the legacy ``utils/somatic_vaf.2.py``. To prove
the strip is output-preserving on the documented mode (no
``--reference`` / ``--conda-env`` passed), this test compares the new
script's output against the golden snapshot committed in M1b.

The test SKIPS until ``tests/data/golden/filtering/somatic_vaf/`` exists,
which happens after M1b runs on the cluster. M2 may land first; the test
activates automatically when the goldens are added.
"""

from __future__ import annotations

from pathlib import Path

import pytest

GOLDEN_DIR = Path(__file__).parent / "data" / "golden" / "filtering" / "somatic_vaf"

pytestmark = pytest.mark.skipif(
    not GOLDEN_DIR.exists(),
    reason="M1b golden snapshots not yet generated — run on cluster",
)


def test_somatic_vaf_filter_matches_golden(tmp_path: Path) -> None:
    """Byte-equivalent (with 1e-6 tolerance on AF columns) against M1b golden.

    Expected fixture layout once M1b lands:
        tests/data/golden/filtering/somatic_vaf/
            input/snv_list.txt        (4-column TSV)
            input/sample.bam          (or symlink into tests/data/)
            expected/output.tsv       (full somatic_vaf.2.py output)
            cmd.txt                   (exact argv used to generate)
    """
    # The detailed assertion logic will be authored when goldens are
    # available so the comparator (numeric-tolerance on the VAF/p-binom
    # columns vs byte-equality on chrom/pos/ref/alt/depth) can be wired
    # against real data. For now the test exists to document the contract
    # and to assert that the harness wakes up once the fixtures land.
    raise AssertionError(
        "M1b goldens present but characterization test body not yet authored"
    )
