"""
test_auto_params.py — Unit tests for scripts/auto_params.py

Coverage targets:
    get_sort_memory()          — fixed 8 GB buffer, clamped to RAM/4
    get_bqsr_memory_gb()       — fixed 8 GB heap, clamped to RAM/4, floor 4 GB
    get_markdup_memory_gb()    — fixed 16 GB heap, clamped to RAM/4, floor 4 GB
    get_gatk_memory_gb()       — fixed 8 GB heap, clamped to RAM/4, floor 4 GB
    get_bwa_threads()          — cpu_count, capped at 4 on low-memory systems
    detect_sequencer()         — instrument ID pattern matching
    get_optical_duplicate_pixel_distance() — patterned vs unpatterned flowcell
"""

from unittest.mock import MagicMock, patch

import pytest
from scripts.auto_params import (
    detect_sequencer,
    get_bqsr_memory_gb,
    get_bwa_threads,
    get_gatk_memory_gb,
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
    m.total = total_gb << 30  # bytes; (>> 30) == total_gb
    m.available = avail_gb << 30  # bytes; (>> 20) == avail_gb * 1024 MB
    return m


def _resolve_with_ram(resolver, total_gb: int):
    """Call ``resolver()`` on a host mocked to have ``total_gb`` of total RAM."""
    with patch(
        "scripts.auto_params.psutil.virtual_memory",
        return_value=_mem(total_gb, total_gb * 3 // 4),
    ):
        return resolver()


# ---------------------------------------------------------------------------
# get_sort_memory
# ---------------------------------------------------------------------------


class TestGetSortMemory:
    """sambamba sort buffer: fixed 8 GB, clamped to a quarter of total RAM."""

    def test_ample_ram_returns_fixed_buffer(self):
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(98, 60)):
            assert get_sort_memory() == "8GB"

    def test_small_host_clamped_to_quarter_of_total(self):
        # 16 GB total → quarter = 4 GB, below the 8 GB default
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(16, 12)):
            assert get_sort_memory() == "4GB"

    def test_tiny_host_floored_at_2gb(self):
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(4, 3)):
            assert get_sort_memory() == "2GB"

    def test_psutil_failure_returns_default(self):
        with patch("scripts.auto_params.psutil.virtual_memory", side_effect=RuntimeError):
            assert get_sort_memory() == "8GB"

    def test_independent_of_available_ram(self):
        """Two hosts with identical total RAM must resolve identically.

        The buffer used to be 30% of *available* RAM, which made the resolved
        parameters depend on whatever else happened to be running at resolve
        time — the same cohort could sort with different buffers on re-run.
        """
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(98, 90)):
            idle_host = get_sort_memory()
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(98, 3)):
            busy_host = get_sort_memory()
        assert idle_host == busy_host


# ---------------------------------------------------------------------------
# get_bqsr_memory_gb / get_gatk_memory_gb / get_markdup_memory_gb
#
# These feed each rule's `resources.mem_mb`, which Snakemake treats as an
# admission-control budget (concurrency = --resources mem_mb // per-job mem_mb).
# They must therefore stay FLAT as total RAM grows: a bigger host has to buy
# more concurrent jobs, not a bigger heap per job.
# ---------------------------------------------------------------------------


class TestGetBqsrMemoryGb:
    """GATK streaming-walker heap: fixed 8 GB, clamped to RAM/4, floor 4 GB."""

    def test_normal_system(self):
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(98, 60)):
            assert get_bqsr_memory_gb() == 8

    def test_low_memory_returns_minimum_4gb(self):
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(4, 3)):
            assert get_bqsr_memory_gb() == 4

    def test_psutil_failure_returns_default(self):
        with patch("scripts.auto_params.psutil.virtual_memory", side_effect=RuntimeError):
            assert get_bqsr_memory_gb() == 8

    def test_does_not_grow_with_total_ram(self):
        """Regression: heap scaled to total RAM starved the scatter rules.

        On the 180-core/98 GB host, `total_gb // 2` produced a 49 GB
        reservation against ~13 GB actual RSS, so `mutect2_scatter` and
        `apply_bqsr` ran ONE job at a time under `--resources mem_mb=90000`
        and 177 of 180 cores sat idle.
        """
        sizes = {gb: _resolve_with_ram(get_bqsr_memory_gb, gb) for gb in (64, 98, 128, 512)}
        assert len(set(sizes.values())) == 1, (
            f"heap must not scale with total RAM, got {sizes}"
        )


class TestGetGatkMemoryGb:
    """General GATK-tool heap: same fixed 8 GB contract as the BQSR heap."""

    def test_normal_system(self):
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(98, 60)):
            assert get_gatk_memory_gb() == 8

    def test_does_not_grow_with_total_ram(self):
        sizes = {gb: _resolve_with_ram(get_gatk_memory_gb, gb) for gb in (32, 98, 512)}
        assert len(set(sizes.values())) == 1, f"got {sizes}"


class TestGetMarkdupMemoryGb:
    """MarkDuplicates heap: fixed 16 GB, clamped to RAM/4, floor 4 GB."""

    def test_normal_system(self):
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(98, 60)):
            assert get_markdup_memory_gb() == 16

    def test_small_host_clamped_to_quarter_of_total(self):
        # 32 GB total → quarter = 8 GB, below the 16 GB default
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(32, 24)):
            assert get_markdup_memory_gb() == 8

    def test_low_memory_returns_minimum_4gb(self):
        with patch("scripts.auto_params.psutil.virtual_memory", return_value=_mem(4, 3)):
            assert get_markdup_memory_gb() == 4

    def test_psutil_failure_returns_default(self):
        with patch("scripts.auto_params.psutil.virtual_memory", side_effect=RuntimeError):
            assert get_markdup_memory_gb() == 16

    def test_does_not_grow_with_total_ram(self):
        # 64 GB and up all clamp-free, so the value must be identical.
        sizes = {gb: _resolve_with_ram(get_markdup_memory_gb, gb) for gb in (64, 98, 512)}
        assert len(set(sizes.values())) == 1, f"got {sizes}"


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

    @pytest.mark.parametrize(
        "instrument,expected",
        [
            ("LH00204", "NovaSeq X"),
            ("A00100", "NovaSeq 6000"),
            ("A00266", "NovaSeq 6000"),
            ("E00143", "HiSeq X"),
            ("K00145", "HiSeq X"),
            ("J00120", "HiSeq 3000/4000"),
            ("SN0196", "HiSeq 2500"),
            ("D00195", "HiSeq 2500"),
            ("HWI-ST123", "HiSeq 2500"),
            ("M03213", "MiSeq"),
            ("NS500487", "NextSeq 500/550"),
            ("NB501234", "NextSeq 500/550"),
            ("VH00204", "NextSeq 2000"),
        ],
    )
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

    @pytest.mark.parametrize(
        "sequencer,expected",
        [
            ("NovaSeq 6000", 2500),
            ("NovaSeq X", 2500),
            ("HiSeq X", 2500),
            ("NovaSeq", 2500),
            ("HiSeq 2500", 100),
            ("HiSeq 3000/4000", 100),
            ("MiSeq", 100),
            ("NextSeq 500/550", 100),
            ("NextSeq 2000", 100),
            ("Unknown", 100),
        ],
    )
    def test_odpd_by_sequencer(self, sequencer, expected):
        assert get_optical_duplicate_pixel_distance(sequencer) == expected
