"""
test_library.py — Unit tests for library/ modules.

Coverage targets:
    library/misc.py      — coroutine decorator, printer
    library/parser.py    — filetype, sample_list
    library/pileup.py    — bases_clean, base_n coroutine, base_qual coroutine
    library/config.py    — run_info_append, log_dir, save_hold_jid
    library/job_queue.py — GridEngineQueue (subprocess mocked)

Excluded (require external processes):
    library/config.py::read_config      — needs `conda info` subprocess
    library/config.py::run_info         — depends on read_config
    library/pileup.py::pileup           — needs samtools subprocess
    library/pileup.py::load_config      — depends on read_config
"""

import os
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

PROJECT_ROOT = Path(__file__).parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


# =============================================================================
# library/misc.py
# =============================================================================

class TestCoroutine:
    """Tests for the @coroutine decorator in library/misc.py."""

    def test_coroutine_is_preadvanced(self):
        """Decorated generator should be at first yield before any send()."""
        from library.misc import coroutine

        @coroutine
        def echo():
            result = None
            while True:
                value = (yield result)
                result = value * 2

        gen = echo()
        assert gen.send(5) == 10

    def test_coroutine_accepts_args(self):
        """Decorator must forward constructor arguments to the generator."""
        from library.misc import coroutine

        @coroutine
        def multiplier(factor):
            result = None
            while True:
                value = (yield result)
                result = value * factor

        gen = multiplier(3)
        assert gen.send(4) == 12

    def test_coroutine_returns_generator_with_send(self):
        """Return value must expose .send() so it can be driven."""
        from library.misc import coroutine

        @coroutine
        def noop():
            while True:
                yield (yield)

        gen = noop()
        assert hasattr(gen, "send")

    def test_coroutine_multiple_sends(self):
        """Each send() should produce an independent result."""
        from library.misc import coroutine

        @coroutine
        def echo():
            result = None
            while True:
                value = (yield result)
                result = value

        gen = echo()
        assert gen.send("first") == "first"
        assert gen.send("second") == "second"


class TestPrinter:
    """Tests for printer() in library/misc.py."""

    def test_prints_to_stdout(self, capsys):
        from library.misc import printer
        printer("hello world")
        captured = capsys.readouterr()
        assert "hello world" in captured.out

    def test_prints_with_flush(self, capsys):
        from library.misc import printer
        printer(42)
        captured = capsys.readouterr()
        assert "42" in captured.out

    def test_broken_pipe_does_not_raise(self):
        """BrokenPipeError on stdout.close() must be swallowed."""
        from library.misc import printer
        with patch("builtins.print", side_effect=BrokenPipeError):
            with patch("sys.stdout") as mock_out, patch("sys.stderr") as mock_err:
                mock_out.close.side_effect = BrokenPipeError
                mock_err.close.side_effect = BrokenPipeError
                # Should not propagate any exception
                printer("test")


# =============================================================================
# library/parser.py
# =============================================================================

class TestFiletype:
    """Tests for filetype() in library/parser.py."""

    @pytest.mark.parametrize("fname,expected", [
        ("sample.bam",     "bam"),
        ("sample.bai",     "bam"),
        ("sample.cram",    "cram"),
        ("sample.crai",    "cram"),
        ("sample.fastq",   "fastq"),
        ("sample.fq",      "fastq"),
        ("sample.fastq.gz","fastq"),
        ("sample.fq.gz",   "fastq"),
        ("sample.bam.gz",  "bam"),
    ])
    def test_recognized_extensions(self, fname, expected):
        from library.parser import filetype
        assert filetype(fname) == expected

    def test_unknown_extension_raises(self):
        from library.parser import filetype
        with pytest.raises(Exception, match="not allowed filetype"):
            filetype("sample.vcf")

    def test_unknown_extension_message_includes_ext(self):
        from library.parser import filetype
        with pytest.raises(Exception, match=r"\.txt"):
            filetype("sample.txt")


class TestSampleList:
    """Tests for sample_list() in library/parser.py."""

    def test_basic_bam_entry(self, tmp_path):
        from library.parser import sample_list
        sfile = tmp_path / "samples.txt"
        sfile.write_text("S1\tsample.bam\t/data\n")
        result = sample_list(str(sfile))
        assert ("S1", "bam") in result
        assert result[("S1", "bam")] == [("sample.bam", "/data")]

    def test_basic_cram_entry(self, tmp_path):
        from library.parser import sample_list
        sfile = tmp_path / "samples.txt"
        sfile.write_text("S2\tsample.cram\t/data\n")
        result = sample_list(str(sfile))
        assert ("S2", "cram") in result

    def test_comment_lines_skipped(self, tmp_path):
        from library.parser import sample_list
        sfile = tmp_path / "samples.txt"
        sfile.write_text("# header comment\nS1\tsample.bam\t/data\n")
        result = sample_list(str(sfile))
        assert len(result) == 1

    def test_multiple_files_same_sample(self, tmp_path):
        from library.parser import sample_list
        sfile = tmp_path / "samples.txt"
        sfile.write_text(
            "S1\tlib1.bam\t/data\n"
            "S1\tlib2.bam\t/data\n"
        )
        result = sample_list(str(sfile))
        assert len(result[("S1", "bam")]) == 2

    def test_fastq_entry(self, tmp_path):
        from library.parser import sample_list
        sfile = tmp_path / "samples.txt"
        sfile.write_text("S3\treads.fastq\t/data\n")
        result = sample_list(str(sfile))
        assert ("S3", "fastq") in result

    def test_empty_file_returns_empty_dict(self, tmp_path):
        from library.parser import sample_list
        sfile = tmp_path / "samples.txt"
        sfile.write_text("")
        result = sample_list(str(sfile))
        assert len(result) == 0

    def test_only_comments_returns_empty_dict(self, tmp_path):
        from library.parser import sample_list
        sfile = tmp_path / "samples.txt"
        sfile.write_text("# comment1\n# comment2\n")
        result = sample_list(str(sfile))
        assert len(result) == 0


# =============================================================================
# library/pileup.py — pure / coroutine functions
# =============================================================================

class TestBasesClean:
    """Tests for bases_clean() in library/pileup.py."""

    def test_plain_bases_unchanged(self):
        from library.pileup import bases_clean
        assert bases_clean("ACGTacgt") == "ACGTacgt"

    def test_ref_match_dots_and_commas_unchanged(self):
        from library.pileup import bases_clean
        assert bases_clean(".,.,") == ".,.,."[:-1]  # ".,." kept

    def test_removes_read_start_marker(self):
        """^X removes ^ and the following quality character."""
        from library.pileup import bases_clean
        result = bases_clean("^!A")
        assert "^" not in result
        assert "!" not in result
        assert "A" in result

    def test_removes_read_end_marker(self):
        from library.pileup import bases_clean
        result = bases_clean("A$T")
        assert "$" not in result
        assert "A" in result
        assert "T" in result

    def test_removes_deletion(self):
        """Deletion notation -2AA removes the flag and the skipped bases."""
        from library.pileup import bases_clean
        result = bases_clean("A-2AAT")
        assert "-" not in result
        # The two A's that are part of the deletion token should be stripped
        assert result == "AT"

    def test_removes_insertion(self):
        """+3ACG removes the insertion token."""
        from library.pileup import bases_clean
        result = bases_clean("A+3ACGT")
        assert "+" not in result
        assert result == "AT"

    def test_empty_string(self):
        from library.pileup import bases_clean
        assert bases_clean("") == ""

    def test_multiple_markers_in_one_string(self):
        from library.pileup import bases_clean
        result = bases_clean("^!A$T")
        assert "^" not in result
        assert "$" not in result
        assert "A" in result
        assert "T" in result


class TestBaseN:
    """Tests for the base_n() coroutine in library/pileup.py."""

    def test_counts_each_base(self):
        from library.pileup import base_n
        gen = base_n()
        result = gen.send(("AACGTacgt", ""))
        assert result["A"] == 2
        assert result["C"] == 1
        assert result["G"] == 1
        assert result["T"] == 1
        assert result["a"] == 1
        assert result["c"] == 1
        assert result["g"] == 1
        assert result["t"] == 1

    def test_empty_bases_all_zero(self):
        from library.pileup import base_n
        gen = base_n()
        result = gen.send(("", ""))
        for base in ("A", "C", "G", "T", "a", "c", "g", "t"):
            assert result[base] == 0

    def test_counts_deletions(self):
        from library.pileup import base_n
        gen = base_n()
        result = gen.send(("**A", ""))
        assert result["dels"] == 2
        assert result["A"] == 1

    def test_multiple_sends_are_independent(self):
        from library.pileup import base_n
        gen = base_n()
        r1 = gen.send(("AAA", ""))
        assert r1["A"] == 3
        r2 = gen.send(("GGG", ""))
        assert r2["G"] == 3
        assert r2["A"] == 0

    def test_has_all_eight_base_keys(self):
        from library.pileup import base_n
        gen = base_n()
        result = gen.send(("A", ""))
        for base in ("A", "C", "G", "T", "a", "c", "g", "t"):
            assert base in result, f"Missing key '{base}'"


class TestBaseQual:
    """Tests for the base_qual() coroutine in library/pileup.py."""

    def test_returns_list_of_base_qual_tuples(self):
        """'!' = ASCII 33 → phred quality 0."""
        from library.pileup import base_qual
        gen = base_qual()
        result = gen.send(("A", "!"))
        assert result == [("A", 0)]

    def test_lowercase_base_uppercased(self):
        from library.pileup import base_qual
        gen = base_qual()
        result = gen.send(("a", "!"))
        assert result[0][0] == "A"

    def test_quality_encoding(self):
        """'I' = ASCII 73 → phred quality 40."""
        from library.pileup import base_qual
        gen = base_qual()
        result = gen.send(("A", "I"))
        assert result[0][1] == 40

    def test_deletion_stars_removed_before_pairing(self):
        """'*' (deletion) should be stripped so it doesn't pair with a qual."""
        from library.pileup import base_qual
        gen = base_qual()
        result = gen.send(("*A", "!!"))
        assert len(result) == 1
        assert result[0][0] == "A"

    def test_empty_bases_returns_empty_list(self):
        from library.pileup import base_qual
        gen = base_qual()
        result = gen.send(("", ""))
        assert result == []

    def test_multiple_bases(self):
        from library.pileup import base_qual
        gen = base_qual()
        result = gen.send(("ACG", "!\"#"))
        assert len(result) == 3
        assert result[0] == ("A", 0)
        assert result[1] == ("C", 1)
        assert result[2] == ("G", 2)


# =============================================================================
# library/config.py — file I/O helpers only
# =============================================================================

class TestRunInfoAppend:
    """Tests for run_info_append() in library/config.py."""

    def test_writes_line_to_file(self, tmp_path):
        from library.config import run_info_append
        fname = str(tmp_path / "run_info.txt")
        run_info_append(fname, "KEY=value")
        assert "KEY=value\n" in Path(fname).read_text()

    def test_appends_multiple_lines(self, tmp_path):
        from library.config import run_info_append
        fname = str(tmp_path / "run_info.txt")
        run_info_append(fname, "LINE1")
        run_info_append(fname, "LINE2")
        lines = Path(fname).read_text().splitlines()
        assert lines == ["LINE1", "LINE2"]

    def test_creates_file_if_not_exists(self, tmp_path):
        from library.config import run_info_append
        fname = str(tmp_path / "new_run_info.txt")
        run_info_append(fname, "DATA")
        assert Path(fname).exists()


class TestLogDir:
    """Tests for log_dir() in library/config.py."""

    def test_creates_logs_subdirectory(self, tmp_path):
        from library.config import log_dir
        sample = str(tmp_path / "SAMPLE1")
        result = log_dir(sample)
        assert Path(result).exists()

    def test_returned_path_ends_with_logs(self, tmp_path):
        from library.config import log_dir
        sample = str(tmp_path / "SAMPLE1")
        result = log_dir(sample)
        assert result == sample + "/logs"

    def test_idempotent_when_called_twice(self, tmp_path):
        from library.config import log_dir
        sample = str(tmp_path / "SAMPLE1")
        result1 = log_dir(sample)
        result2 = log_dir(sample)
        assert result1 == result2
        assert Path(result1).exists()


class TestSaveHoldJid:
    """Tests for save_hold_jid() in library/config.py."""

    def test_writes_jid_to_file(self, tmp_path):
        from library.config import save_hold_jid
        fname = str(tmp_path / "run.jid")
        save_hold_jid(fname, "12345")
        content = Path(fname).read_text().strip()
        assert content == "12345"

    def test_creates_parent_directories(self, tmp_path):
        from library.config import save_hold_jid
        fname = str(tmp_path / "deep" / "nested" / "run.jid")
        save_hold_jid(fname, "99999")
        assert Path(fname).exists()

    def test_overwrites_existing_file(self, tmp_path):
        from library.config import save_hold_jid
        fname = str(tmp_path / "run.jid")
        save_hold_jid(fname, "111")
        save_hold_jid(fname, "222")
        content = Path(fname).read_text().strip()
        assert content == "222"


# =============================================================================
# library/job_queue.py — subprocess mocked
# =============================================================================

class TestGridEngineQueueInit:

    def test_run_jid_initially_none(self):
        from library.job_queue import GridEngineQueue
        q = GridEngineQueue()
        assert q.run_jid is None


class TestSetRunJid:

    def test_stores_path(self, tmp_path):
        from library.job_queue import GridEngineQueue
        q = GridEngineQueue()
        fname = str(tmp_path / "run.jid")
        Path(fname).touch()
        q.set_run_jid(fname)
        assert q.run_jid == fname

    def test_new_flag_creates_empty_file(self, tmp_path):
        from library.job_queue import GridEngineQueue
        q = GridEngineQueue()
        subdir = tmp_path / "subdir"
        subdir.mkdir()
        fname = str(subdir / "run.jid")
        q.set_run_jid(fname, new=True)
        assert Path(fname).exists()
        assert Path(fname).read_text() == ""

    def test_new_flag_sets_run_jid(self, tmp_path):
        from library.job_queue import GridEngineQueue
        q = GridEngineQueue()
        fname = str(tmp_path / "run.jid")
        q.set_run_jid(fname, new=True)
        assert q.run_jid == fname


class TestNumRunJidInQueue:

    def test_returns_zero_if_file_missing(self):
        from library.job_queue import GridEngineQueue
        q = GridEngineQueue()
        assert q.num_run_jid_in_queue("/nonexistent/path.jid") == 0

    def test_returns_zero_if_file_empty(self, tmp_path):
        from library.job_queue import GridEngineQueue
        q = GridEngineQueue()
        fname = tmp_path / "run.jid"
        fname.write_text("")
        assert q.num_run_jid_in_queue(str(fname)) == 0

    def test_queries_squeue_for_active_jobs(self, tmp_path):
        from library.job_queue import GridEngineQueue
        q = GridEngineQueue()
        fname = tmp_path / "run.jid"
        fname.write_text("12345\n")
        mock_result = MagicMock()
        mock_result.stdout = "1\n"
        with patch("subprocess.run", return_value=mock_result):
            count = q.num_run_jid_in_queue(str(fname))
        assert count == 1

    def test_returns_zero_when_squeue_reports_none(self, tmp_path):
        from library.job_queue import GridEngineQueue
        q = GridEngineQueue()
        fname = tmp_path / "run.jid"
        fname.write_text("12345\n")
        mock_result = MagicMock()
        mock_result.stdout = "0\n"
        with patch("subprocess.run", return_value=mock_result):
            count = q.num_run_jid_in_queue(str(fname))
        assert count == 0


class TestSubmit:

    def test_submit_calls_sbatch(self):
        from library.job_queue import GridEngineQueue
        q = GridEngineQueue()
        mock_result = MagicMock()
        mock_result.stdout = "Submitted batch job 42"
        with patch("subprocess.run", return_value=mock_result) as mock_run:
            q.submit("-p queue", "my_script.sh")
        cmd = mock_run.call_args[0][0]
        assert cmd[0] == "sbatch"

    def test_submit_returns_stdout(self):
        from library.job_queue import GridEngineQueue
        q = GridEngineQueue()
        mock_result = MagicMock()
        mock_result.stdout = "Submitted batch job 42"
        with patch("subprocess.run", return_value=mock_result):
            jid = q.submit("-p queue", "my_script.sh")
        assert jid == "Submitted batch job 42"

    def test_submit_appends_jid_to_run_jid_file(self, tmp_path):
        from library.job_queue import GridEngineQueue
        q = GridEngineQueue()
        fname = tmp_path / "run.jid"
        fname.touch()
        q.set_run_jid(str(fname))
        mock_result = MagicMock()
        mock_result.stdout = "42"
        with patch("subprocess.run", return_value=mock_result):
            q.submit("", "script.sh")
        assert "42" in fname.read_text()

    def test_submit_without_run_jid_does_not_fail(self):
        """run_jid=None means no file to append to; submit should still work."""
        from library.job_queue import GridEngineQueue
        q = GridEngineQueue()
        mock_result = MagicMock()
        mock_result.stdout = "99"
        with patch("subprocess.run", return_value=mock_result):
            jid = q.submit("", "script.sh")
        assert jid == "99"
