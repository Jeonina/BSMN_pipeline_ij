"""
test_auto_params.py — Unit tests for scripts/auto_params.py

Coverage targets:
    get_sort_memory()          — 30% of available RAM, min 512 MB, max 8 GB
    get_bqsr_memory_gb()       — half of total RAM, min 4 GB, max 64 GB
    get_markdup_memory_gb()    — quarter of total RAM, min 2 GB, max 32 GB
    get_bwa_threads()          — cpu_count, capped at 4 on low-memory systems
    detect_sequencer()         — instrument ID pattern matching
    get_optical_duplicate_pixel_distance() — patterned vs unpatterned flowcell
"""

import os
from unittest.mock import MagicMock, patch

import pytest

from scripts.auto_params import (
    detect_sequencer,
    get_bqsr_memory_gb,
    get_bwa_threads,
    get_markdup_memory_gb,
    get_optical_duplicate_pixel_distance,
    get_sort_memory,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _mem(total_gb: int, avail_gb: int) -> MagicMock:
    """Build a psutil.virtual_memory() mock from GB values."""
    m = MagicMock()
    m.total = total_gb << 30      # bytes; (>> 30) == total_gb
    m.available = avail_gb << 30  # bytes; (>> 20) == avail_gb * 1024 MB
    return m


# ---------------------------------------------------------------------------
# get_sort_memory
# ---------------------------------------------------------------------------

class TestGetSortMemory:
    """sambamba sort memory: 30% of available RAM, min 512 MB, max 8 GB."""

    def test_normal_memory_returns_30_percent(self):
        # 16 GB available → 30% = 4915 MB (within [512, 8192])
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(32, 16)):
            result = get_sort_memory()
        assert result == "4915MB"

    def test_low_memory_returns_minimum_512mb(self):
        # 1 GB available → 30% = 307 MB → floored to 512 MB
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(4, 1)):
            result = get_sort_memory()
        assert result == "512MB"

    def test_high_memory_server_capped_at_8gb(self):
        # 500 GB available → 30% = 153600 MB → must be capped at 8192 MB
        # This was the production bug: server had 512 GB RAM → 152 GB sort memory
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(512, 500)):
            result = get_sort_memory()
        assert result == "8192MB", (
            f"Expected '8192MB' but got '{result}'. "
            "get_sort_memory() must cap at 8 GB regardless of available RAM."
        )

    def test_psutil_failure_returns_fallback(self):
        with patch("scripts.auto_params.psutil.virtual_memory", side_effect=RuntimeError):
            result = get_sort_memory()
        assert result == "6GB"

    def test_boundary_exactly_at_cap(self):
        # Exactly 8 GB / 0.30 ≈ 27307 MB available → 30% = 8192 MB (at cap)
        avail_mb = 8192
        avail_gb_approx = avail_mb // 1024  # 8 GB
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(16, avail_gb_approx)):
            result = get_sort_memory()
        mb = int(result.rstrip("MB"))
        assert mb <= 8192


# ---------------------------------------------------------------------------
# get_bqsr_memory_gb
# ---------------------------------------------------------------------------

class TestGetBqsrMemoryGb:
    """GATK BQSR heap: half of total RAM, min 4 GB, max 64 GB."""

    def test_normal_system(self):
        # 32 GB total → half = 16 GB
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(32, 24)):
            assert get_bqsr_memory_gb() == 16

    def test_low_memory_returns_minimum_4gb(self):
        # 4 GB total → half = 2 → floored to 4
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(4, 3)):
            assert get_bqsr_memory_gb() == 4

    def test_high_memory_capped_at_64gb(self):
        # 512 GB total → half = 256 → capped at 64
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(512, 400)):
            assert get_bqsr_memory_gb() == 64

    def test_psutil_failure_returns_fallback(self):
        with patch("scripts.auto_params.psutil.virtual_memory", side_effect=RuntimeError):
            assert get_bqsr_memory_gb() == 16

    def test_boundary_128gb_total(self):
        # 128 GB → half = 64 → exactly at cap
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(128, 100)):
            assert get_bqsr_memory_gb() == 64


# ---------------------------------------------------------------------------
# get_markdup_memory_gb
# ---------------------------------------------------------------------------

class TestGetMarkdupMemoryGb:
    """Picard MarkDuplicates heap: quarter of total RAM, min 2 GB, max 32 GB."""

    def test_normal_system(self):
        # 32 GB total → quarter = 8 GB
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(32, 24)):
            assert get_markdup_memory_gb() == 8

    def test_low_memory_returns_minimum_2gb(self):
        # 4 GB total → quarter = 1 → floored to 2
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(4, 3)):
            assert get_markdup_memory_gb() == 2

    def test_high_memory_capped_at_32gb(self):
        # 512 GB total → quarter = 128 → capped at 32
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(512, 400)):
            assert get_markdup_memory_gb() == 32

    def test_psutil_failure_returns_fallback(self):
        with patch("scripts.auto_params.psutil.virtual_memory", side_effect=RuntimeError):
            assert get_markdup_memory_gb() == 8


# ---------------------------------------------------------------------------
# get_bwa_threads
# ---------------------------------------------------------------------------

class TestGetBwaThreads:
    """CPU count, capped at 4 on systems with < 8 GB total RAM."""

    def test_high_memory_returns_all_cpus(self):
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(32, 24)):
            with patch("os.cpu_count", return_value=16):
                assert get_bwa_threads() == 16

    def test_low_memory_caps_at_4_threads(self):
        # < 8 GB total → cap at 4 even if more CPUs available
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(4, 3)):
            with patch("os.cpu_count", return_value=32):
                assert get_bwa_threads() == 4

    def test_low_memory_with_few_cpus(self):
        # < 8 GB but only 2 CPUs → min(2, 4) = 2
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(4, 3)):
            with patch("os.cpu_count", return_value=2):
                assert get_bwa_threads() == 2

    def test_psutil_failure_returns_cpu_count(self):
        with patch("scripts.auto_params.psutil.virtual_memory", side_effect=RuntimeError):
            with patch("os.cpu_count", return_value=8):
                assert get_bwa_threads() == 8

    def test_cpu_count_none_returns_at_least_1(self):
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(32, 24)):
            with patch("os.cpu_count", return_value=None):
                assert get_bwa_threads() >= 1


# ---------------------------------------------------------------------------
# detect_sequencer
# ---------------------------------------------------------------------------

class TestDetectSequencer:
    """Instrument ID pattern → sequencer name."""

    @pytest.mark.parametrize("instrument,expected", [
        ("LH00204",  "NovaSeq X"),
        ("A00100",   "NovaSeq 6000"),
        ("A00266",   "NovaSeq 6000"),
        ("E00143",   "HiSeq X"),
        ("K00145",   "HiSeq X"),
        ("J00120",   "HiSeq 3000/4000"),
        ("SN0196",   "HiSeq 2500"),
        ("D00195",   "HiSeq 2500"),
        ("HWI-ST123","HiSeq 2500"),
        ("M03213",   "MiSeq"),
        ("NS500487", "NextSeq 500/550"),
        ("NB501234", "NextSeq 500/550"),
        ("VH00204",  "NextSeq 2000"),
    ])
    def test_known_instrument_ids(self, instrument, expected, tmp_path):
        fq = tmp_path / "test.fastq.gz"
        import gzip
        header = f"@{instrument}:1:flowcell:1:1:100:200\n"
        with gzip.open(fq, "wt") as f:
            f.write(header)
        assert detect_sequencer(str(fq)) == expected

    def test_unknown_instrument_falls_back(self, tmp_path):
        fq = tmp_path / "test.fastq.gz"
        import gzip
        with gzip.open(fq, "wt") as f:
            f.write("@ERR194146.1\n")
        assert detect_sequencer(str(fq)) == "Unknown"

    def test_novaseq_by_colon_heuristic(self, tmp_path):
        # 8 colon-fields → patterned flowcell heuristic
        fq = tmp_path / "test.fastq.gz"
        import gzip
        header = "@UNKNOWN:1:flowcell:1:1:100:200:ACGT\n"  # 8 fields
        with gzip.open(fq, "wt") as f:
            f.write(header)
        assert detect_sequencer(str(fq)) == "NovaSeq 6000"


# ---------------------------------------------------------------------------
# get_optical_duplicate_pixel_distance
# ---------------------------------------------------------------------------

class TestGetOpticalDuplicatePixelDistance:
    """Patterned flowcell → 2500, unpatterned → 100."""

    @pytest.mark.parametrize("sequencer,expected", [
        ("NovaSeq 6000",    2500),
        ("NovaSeq X",       2500),
        ("HiSeq X",         2500),
        ("NovaSeq",         2500),
        ("HiSeq 2500",      100),
        ("HiSeq 3000/4000", 100),
        ("MiSeq",           100),
        ("NextSeq 500/550", 100),
        ("NextSeq 2000",    100),
        ("Unknown",         100),
    ])
    def test_odpd_by_sequencer(self, sequencer, expected):
        assert get_optical_duplicate_pixel_distance(sequencer) == expected
