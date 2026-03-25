"""
TDD test suite for BSMN mapping pipeline (Phase 1).

Test order follows PIPELINE.md spec:
  1. auto_params unit tests   — sequencer detection, ODPD, system resources
  2. containers.yaml validation
  3. Snakemake dry-run (DAG check)  — marked @pytest.mark.integration

Run all:
    pytest tests/test_mapping.py -v

Skip integration tests (no snakemake install required):
    pytest tests/test_mapping.py -v -m "not integration"
"""

import gzip
import os
import subprocess
import tempfile
from pathlib import Path

import pytest
import yaml

PROJECT_ROOT = Path(__file__).parent.parent
WORKFLOW_DIR = PROJECT_ROOT / "workflow"
CONFIG_DIR = PROJECT_ROOT / "config"
TESTS_DIR = PROJECT_ROOT / "tests"
DATA_DIR = TESTS_DIR / "data"


# =============================================================================
# Helpers
# =============================================================================

def _make_fastq(path: Path, instrument: str, n_fields: int, flowcell: str = "AABBCCDD") -> str:
    """
    Write a tiny synthetic FASTQ with the given instrument ID and field count.

    n_fields: total colon-separated fields in the read name (e.g. 8 → 7 colons).
    CASAVA 1.8+ layout: instrument:run:flowcell:lane:tile:x:y[(:umi)]
    """
    extra_fields = ["1"] * max(0, n_fields - 7)  # pad beyond standard 7
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as fh:
        for i in range(20):
            # Build read name with exactly n_fields colon-separated parts
            parts = [instrument, "123", flowcell, "1", "1101",
                     str(1000 + i), str(2000 + i)] + extra_fields
            parts = parts[:n_fields]
            rname = ":".join(parts)
            fh.write(f"@{rname} 1:N:0:ATCG\n")
            fh.write("ACGTACGTACGTACGTACGTACGTACGTACGT\n")
            fh.write("+\n")
            fh.write("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")
    return str(path)


# =============================================================================
# Session fixtures — synthetic FASTQ per sequencer model
# =============================================================================

@pytest.fixture(scope="session")
def novaseq_fastq(tmp_path_factory):
    """NovaSeq 6000: instrument A00100, 8 fields (7 colons)."""
    p = tmp_path_factory.mktemp("fastq") / "novaseq6000.R1.fastq.gz"
    return _make_fastq(p, "A00100", n_fields=8)


@pytest.fixture(scope="session")
def novaseq_x_fastq(tmp_path_factory):
    """NovaSeq X: instrument LH00204, 8 fields."""
    p = tmp_path_factory.mktemp("fastq") / "novaseq_x.R1.fastq.gz"
    return _make_fastq(p, "LH00204", n_fields=8)


@pytest.fixture(scope="session")
def hiseq_x_fastq(tmp_path_factory):
    """HiSeq X: instrument E00143, 7 fields (6 colons)."""
    p = tmp_path_factory.mktemp("fastq") / "hiseq_x.R1.fastq.gz"
    return _make_fastq(p, "E00143", n_fields=7)


@pytest.fixture(scope="session")
def hiseq_34k_fastq(tmp_path_factory):
    """HiSeq 3000/4000: instrument J00120, 7 fields."""
    p = tmp_path_factory.mktemp("fastq") / "hiseq_34k.R1.fastq.gz"
    return _make_fastq(p, "J00120", n_fields=7)


@pytest.fixture(scope="session")
def hiseq_fastq(tmp_path_factory):
    """HiSeq 2500: instrument SN0196, 7 fields."""
    p = tmp_path_factory.mktemp("fastq") / "hiseq2500.R1.fastq.gz"
    return _make_fastq(p, "SN0196", n_fields=7)


@pytest.fixture(scope="session")
def miseq_fastq(tmp_path_factory):
    """MiSeq: instrument M03213, 7 fields."""
    p = tmp_path_factory.mktemp("fastq") / "miseq.R1.fastq.gz"
    return _make_fastq(p, "M03213", n_fields=7)


@pytest.fixture(scope="session")
def nextseq_fastq(tmp_path_factory):
    """NextSeq 500/550: instrument NS500487, 7 fields."""
    p = tmp_path_factory.mktemp("fastq") / "nextseq.R1.fastq.gz"
    return _make_fastq(p, "NS500487", n_fields=7)


@pytest.fixture(scope="session")
def nextseq2k_fastq(tmp_path_factory):
    """NextSeq 2000: instrument VH00204, 7 fields."""
    p = tmp_path_factory.mktemp("fastq") / "nextseq2k.R1.fastq.gz"
    return _make_fastq(p, "VH00204", n_fields=7)


# =============================================================================
# 1. auto_params — sequencer detection
# =============================================================================

class TestDetectSequencer:

    # --- patterned flowcell (ODPD=2500) ---

    def test_novaseq_6000_instrument_id(self, novaseq_fastq):
        from scripts.auto_params import detect_sequencer
        assert detect_sequencer(novaseq_fastq) == "NovaSeq 6000"

    def test_novaseq_x_instrument_id(self, novaseq_x_fastq):
        from scripts.auto_params import detect_sequencer
        assert detect_sequencer(novaseq_x_fastq) == "NovaSeq X"

    def test_hiseq_x_instrument_id(self, hiseq_x_fastq):
        from scripts.auto_params import detect_sequencer
        assert detect_sequencer(hiseq_x_fastq) == "HiSeq X"

    # --- unpatterned flowcell (ODPD=100) ---

    def test_hiseq_3000_4000_instrument_id(self, hiseq_34k_fastq):
        from scripts.auto_params import detect_sequencer
        assert detect_sequencer(hiseq_34k_fastq) == "HiSeq 3000/4000"

    def test_hiseq_2500_instrument_id(self, hiseq_fastq):
        from scripts.auto_params import detect_sequencer
        assert detect_sequencer(hiseq_fastq) == "HiSeq 2500"

    def test_miseq_instrument_id(self, miseq_fastq):
        from scripts.auto_params import detect_sequencer
        assert detect_sequencer(miseq_fastq) == "MiSeq"

    def test_nextseq_500_550_instrument_id(self, nextseq_fastq):
        from scripts.auto_params import detect_sequencer
        assert detect_sequencer(nextseq_fastq) == "NextSeq 500/550"

    def test_nextseq_2000_instrument_id(self, nextseq2k_fastq):
        from scripts.auto_params import detect_sequencer
        assert detect_sequencer(nextseq2k_fastq) == "NextSeq 2000"

    # --- fallback behavior ---

    def test_uncompressed_fastq_works(self, tmp_path):
        """detect_sequencer handles plain (non-gz) FASTQ."""
        from scripts.auto_params import detect_sequencer
        fq = tmp_path / "plain.fastq"
        fq.write_text(
            "@A00100:123:AABBCCDD:1:1101:1000:2000:NNNNN 1:N:0:ATCG\n"
            "ACGT\n+\nIIII\n"
        )
        assert detect_sequencer(str(fq)) == "NovaSeq 6000"

    def test_boundary_8_fields_unknown_instrument_is_novaseq(self, tmp_path):
        """8 colon-fields (7 colons) with unrecognized instrument → colon heuristic → NovaSeq 6000."""
        from scripts.auto_params import detect_sequencer
        fq = tmp_path / "boundary.fastq.gz"
        # 8 fields = 7 colons: UNKWN:b:c:d:e:f:g:h
        with gzip.open(fq, "wt") as f:
            f.write("@UNKWN:b:c:d:e:f:g:h 1:N:0:\nACGT\n+\nIIII\n")
        assert detect_sequencer(str(fq)) == "NovaSeq 6000"

    def test_boundary_7_fields_unknown_instrument_is_unknown(self, tmp_path):
        """7 colon-fields (6 colons) with unrecognized instrument → unknown fallback."""
        from scripts.auto_params import detect_sequencer
        fq = tmp_path / "hiseq6.fastq.gz"
        with gzip.open(fq, "wt") as f:
            f.write("@UNKWN:b:c:d:e:f:g 1:N:0:\nACGT\n+\nIIII\n")
        assert detect_sequencer(str(fq)) == "Unknown"

    def test_unknown_instrument_defaults_to_unknown(self, tmp_path):
        """Completely unrecognized instrument ID → Unknown."""
        from scripts.auto_params import detect_sequencer
        fq = tmp_path / "unknown.fastq.gz"
        with gzip.open(fq, "wt") as f:
            f.write("@CUSTOMSEQ:1:FC:1:1101:100:200 1:N:0:\nACGT\n+\nIIII\n")
        assert detect_sequencer(str(fq)) == "Unknown"

    # --- NB prefix (NextSeq alternative naming) ---

    def test_nb_prefix_is_nextseq(self, tmp_path):
        """NB-prefixed instrument (NextSeq 550) → NextSeq 500/550."""
        from scripts.auto_params import detect_sequencer
        fq = tmp_path / "nb.fastq.gz"
        with gzip.open(fq, "wt") as f:
            f.write("@NB501234:123:AABBCCDD:1:1101:1000:2000 1:N:0:\nACGT\n+\nIIII\n")
        assert detect_sequencer(str(fq)) == "NextSeq 500/550"

    # --- HWI prefix (legacy HiSeq 2500) ---

    def test_hwi_prefix_is_hiseq_2500(self, tmp_path):
        """HWI-prefixed instrument (legacy HiSeq) → HiSeq 2500."""
        from scripts.auto_params import detect_sequencer
        fq = tmp_path / "hwi.fastq.gz"
        with gzip.open(fq, "wt") as f:
            f.write("@HWI-ST1276:71:D1B67ACXX:7:1101:1428:89553 1:N:0:\nACGT\n+\nIIII\n")
        assert detect_sequencer(str(fq)) == "HiSeq 2500"


# =============================================================================
# 2. auto_params — sequencer_evidence output
# =============================================================================

class TestSequencerEvidence:

    def test_evidence_keys_present(self, novaseq_fastq, tmp_path):
        """resolve_params must include sequencer_evidence with required keys."""
        from scripts.auto_params import resolve_params
        params = resolve_params(novaseq_fastq, str(tmp_path / "r.yaml"))
        ev = params["sequencer_evidence"]
        for key in ("instrument_id", "colon_fields", "flowcell_id",
                    "detection_method", "read_name_example"):
            assert key in ev, f"sequencer_evidence missing key: {key}"

    def test_novaseq_evidence_method_is_instrument_id(self, novaseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        params = resolve_params(novaseq_fastq, str(tmp_path / "r.yaml"))
        assert params["sequencer_evidence"]["detection_method"] == "instrument_id_pattern"

    def test_novaseq_evidence_instrument_id(self, novaseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        params = resolve_params(novaseq_fastq, str(tmp_path / "r.yaml"))
        assert params["sequencer_evidence"]["instrument_id"] == "A00100"

    def test_colon_fallback_evidence_method(self, tmp_path):
        """8-field unknown instrument → colon_count_heuristic in evidence."""
        from scripts.auto_params import resolve_params
        fq = tmp_path / "colon.fastq.gz"
        with gzip.open(fq, "wt") as f:
            f.write("@UNKWN:b:c:d:e:f:g:h 1:N:0:\nACGT\n+\nIIII\n")
        params = resolve_params(str(fq), str(tmp_path / "r.yaml"))
        assert params["sequencer_evidence"]["detection_method"] == "colon_count_heuristic"

    def test_unknown_fallback_evidence_method(self, tmp_path):
        """Unrecognized instrument with 6 colons → unknown_fallback in evidence."""
        from scripts.auto_params import resolve_params
        fq = tmp_path / "unk.fastq.gz"
        with gzip.open(fq, "wt") as f:
            f.write("@UNKWN:b:c:d:e:f:g 1:N:0:\nACGT\n+\nIIII\n")
        params = resolve_params(str(fq), str(tmp_path / "r.yaml"))
        assert params["sequencer_evidence"]["detection_method"] == "unknown_fallback"

    def test_evidence_persisted_in_yaml(self, novaseq_fastq, tmp_path):
        """sequencer_evidence must be written to resolved_params.yaml."""
        from scripts.auto_params import resolve_params
        out = str(tmp_path / "r.yaml")
        resolve_params(novaseq_fastq, out)
        with open(out) as fh:
            data = yaml.safe_load(fh)
        assert "sequencer_evidence" in data
        assert isinstance(data["sequencer_evidence"], dict)


# =============================================================================
# 3. auto_params — ODPD selection
# =============================================================================

class TestGetODPD:

    # patterned flowcell → 2500
    def test_novaseq_odpd_is_2500(self):
        from scripts.auto_params import get_optical_duplicate_pixel_distance
        assert get_optical_duplicate_pixel_distance("NovaSeq") == 2500  # legacy compat

    def test_novaseq_6000_odpd_is_2500(self):
        from scripts.auto_params import get_optical_duplicate_pixel_distance
        assert get_optical_duplicate_pixel_distance("NovaSeq 6000") == 2500

    def test_novaseq_x_odpd_is_2500(self):
        from scripts.auto_params import get_optical_duplicate_pixel_distance
        assert get_optical_duplicate_pixel_distance("NovaSeq X") == 2500

    def test_hiseq_x_odpd_is_2500(self):
        from scripts.auto_params import get_optical_duplicate_pixel_distance
        assert get_optical_duplicate_pixel_distance("HiSeq X") == 2500

    # unpatterned flowcell → 100
    def test_hiseq_odpd_is_100(self):
        from scripts.auto_params import get_optical_duplicate_pixel_distance
        assert get_optical_duplicate_pixel_distance("HiSeq") == 100  # legacy compat

    def test_hiseq_2500_odpd_is_100(self):
        from scripts.auto_params import get_optical_duplicate_pixel_distance
        assert get_optical_duplicate_pixel_distance("HiSeq 2500") == 100

    def test_hiseq_3000_4000_odpd_is_100(self):
        from scripts.auto_params import get_optical_duplicate_pixel_distance
        assert get_optical_duplicate_pixel_distance("HiSeq 3000/4000") == 100

    def test_miseq_odpd_is_100(self):
        from scripts.auto_params import get_optical_duplicate_pixel_distance
        assert get_optical_duplicate_pixel_distance("MiSeq") == 100

    def test_nextseq_odpd_is_100(self):
        from scripts.auto_params import get_optical_duplicate_pixel_distance
        assert get_optical_duplicate_pixel_distance("NextSeq 500/550") == 100

    def test_nextseq_2000_odpd_is_100(self):
        from scripts.auto_params import get_optical_duplicate_pixel_distance
        assert get_optical_duplicate_pixel_distance("NextSeq 2000") == 100

    def test_unknown_odpd_is_100(self):
        from scripts.auto_params import get_optical_duplicate_pixel_distance
        assert get_optical_duplicate_pixel_distance("Unknown") == 100


# =============================================================================
# 4. auto_params — system resource queries
# =============================================================================

class TestGetBwaThreads:
    def test_returns_positive_int(self):
        from scripts.auto_params import get_bwa_threads
        t = get_bwa_threads()
        assert isinstance(t, int)
        assert t >= 1


class TestGetBqsrMemory:
    def test_returns_positive_int(self):
        from scripts.auto_params import get_bqsr_memory_gb
        m = get_bqsr_memory_gb()
        assert isinstance(m, int)
        assert m >= 4

    def test_bounded_between_4_and_64(self):
        from scripts.auto_params import get_bqsr_memory_gb
        m = get_bqsr_memory_gb()
        assert 4 <= m <= 64


# =============================================================================
# 5. auto_params — resolve_params() end-to-end
# =============================================================================

class TestResolveParams:
    def test_writes_yaml_file(self, novaseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        out = str(tmp_path / "resolved.yaml")
        resolve_params(novaseq_fastq, out)
        assert os.path.exists(out)

    def test_yaml_is_parseable(self, novaseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        out = str(tmp_path / "resolved.yaml")
        resolve_params(novaseq_fastq, out)
        with open(out) as fh:
            data = yaml.safe_load(fh)
        assert isinstance(data, dict)

    def test_required_keys_present(self, novaseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        params = resolve_params(novaseq_fastq, str(tmp_path / "r.yaml"))
        required = [
            "resolved_at",
            "input_fastq",
            "sequencer",
            "sequencer_evidence",
            "bwa_threads",
            "sort_threads",
            "sort_memory",
            "optical_duplicate_pixel_distance",
            "bqsr_memory_gb",
            "markdup_memory",
        ]
        for key in required:
            assert key in params, f"Missing key in resolved_params: {key}"

    def test_novaseq_sets_odpd_2500(self, novaseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        params = resolve_params(novaseq_fastq, str(tmp_path / "r.yaml"))
        assert params["optical_duplicate_pixel_distance"] == 2500

    def test_hiseq_sets_odpd_100(self, hiseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        params = resolve_params(hiseq_fastq, str(tmp_path / "r.yaml"))
        assert params["optical_duplicate_pixel_distance"] == 100

    def test_miseq_sets_odpd_100(self, miseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        params = resolve_params(miseq_fastq, str(tmp_path / "r.yaml"))
        assert params["optical_duplicate_pixel_distance"] == 100

    def test_nextseq_sets_odpd_100(self, nextseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        params = resolve_params(nextseq_fastq, str(tmp_path / "r.yaml"))
        assert params["optical_duplicate_pixel_distance"] == 100

    def test_creates_parent_directory(self, novaseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        out = str(tmp_path / "nested" / "subdir" / "resolved.yaml")
        resolve_params(novaseq_fastq, out)
        assert os.path.exists(out)

    def test_input_fastq_stored_as_absolute_path(self, novaseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        params = resolve_params(novaseq_fastq, str(tmp_path / "r.yaml"))
        assert os.path.isabs(params["input_fastq"])

    def test_resolved_at_is_iso_format(self, novaseq_fastq, tmp_path):
        from scripts.auto_params import resolve_params
        import datetime
        params = resolve_params(novaseq_fastq, str(tmp_path / "r.yaml"))
        datetime.datetime.fromisoformat(params["resolved_at"])


# =============================================================================
# 6. containers.yaml structure validation
# =============================================================================

class TestContainersYaml:
    """Validate config/containers.yaml before any Snakemake run."""

    REQUIRED_TOOLS = ["bwa", "sambamba", "picard", "gatk", "samtools"]
    REQUIRED_FIELDS = ["name", "version", "uri", "sif"]

    @pytest.fixture
    def containers(self):
        path = CONFIG_DIR / "containers.yaml"
        assert path.exists(), f"config/containers.yaml not found at {path}"
        with open(path) as fh:
            return yaml.safe_load(fh)

    def test_required_tools_present(self, containers):
        for tool in self.REQUIRED_TOOLS:
            assert tool in containers, f"Missing tool in containers.yaml: {tool}"

    def test_each_tool_has_required_fields(self, containers):
        for tool in self.REQUIRED_TOOLS:
            for field in self.REQUIRED_FIELDS:
                assert field in containers[tool], (
                    f"containers.yaml[{tool}] missing field: '{field}'"
                )

    def test_sif_paths_under_containers_dir(self, containers):
        for tool, spec in containers.items():
            assert spec["sif"].startswith("containers/"), (
                f"{tool}.sif should start with 'containers/', got: {spec['sif']}"
            )

    def test_uri_starts_with_docker_or_oras(self, containers):
        for tool, spec in containers.items():
            assert spec["uri"].startswith(("docker://", "oras://", "library://")), (
                f"{tool}.uri has unexpected scheme: {spec['uri']}"
            )


# =============================================================================
# 7. Snakemake DAG dry-run  (integration — requires snakemake installed)
# =============================================================================

@pytest.mark.integration
class TestSnakemakeDryRun:
    """
    Verify the Snakemake DAG can be built without errors.
    Provides a temporary environment with minimal config + test samples.
    """

    @pytest.fixture(scope="class")
    def dry_run_env(self, tmp_path_factory, novaseq_fastq):
        """Create a temporary working directory with all required config files."""
        env = tmp_path_factory.mktemp("snakemake_env")
        cfg = env / "config"
        cfg.mkdir()

        # --- samples.tsv ---
        (cfg / "samples.tsv").write_text(
            "sample_id\treadgroup\tfq1\tfq2\n"
            f"TEST001\tRG1\t{novaseq_fastq}\t{novaseq_fastq}\n"
        )

        # --- containers.yaml (copy from project) ---
        import shutil
        shutil.copy(CONFIG_DIR / "containers.yaml", cfg / "containers.yaml")

        # --- resolved_params.yaml (generated inline) ---
        from scripts.auto_params import resolve_params
        resolve_params(novaseq_fastq, str(cfg / "resolved_params.yaml"))

        # --- minimal config.yaml ---
        dummy_ref = str(env / "ref.fasta")
        (env / "ref.fasta").touch()
        (cfg / "config.yaml").write_text(
            f"samples: config/samples.tsv\n"
            f"ref:\n"
            f"  build: hg38\n"
            f"  fasta: {dummy_ref}\n"
            f"  dict: {dummy_ref.replace('.fasta', '.dict')}\n"
            f"known_sites:\n"
            f"  dbsnp: dummy.vcf.gz\n"
            f"  mills: dummy.vcf.gz\n"
            f"  indels: dummy.vcf.gz\n"
            f"mapping:\n"
            f"  bwa_threads: 2\n"
            f"  sort_threads: 1\n"
            f"  sort_memory: 256M\n"
            f"  markdup_memory: 1G\n"
            f"  output_format: cram\n"
        )
        return env

    def test_dry_run_exits_zero(self, dry_run_env):
        result = subprocess.run(
            [
                "snakemake",
                "--dry-run",
                "--snakefile", str(WORKFLOW_DIR / "Snakefile"),
                "--directory", str(dry_run_env),
                "--cores", "1",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, (
            f"Snakemake dry-run failed.\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    def test_all_mapping_rules_present_in_dag(self, dry_run_env):
        result = subprocess.run(
            [
                "snakemake",
                "--dry-run",
                "--snakefile", str(WORKFLOW_DIR / "Snakefile"),
                "--directory", str(dry_run_env),
                "--cores", "1",
            ],
            capture_output=True,
            text=True,
        )
        combined = result.stdout + result.stderr
        expected_rules = [
            "bwa_mem_sort",
            "merge_bams",
            "mark_duplicates",
            "base_recalibrator",
            "apply_bqsr",
            "samtools_flagstat",
        ]
        for rule in expected_rules:
            assert rule in combined, (
                f"Rule '{rule}' not found in dry-run output.\n"
                f"Output:\n{combined[:2000]}"
            )
