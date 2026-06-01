"""Tests for BAM/CRAM ingest mode (external-alignment input).

Feature: the pipeline accepts a pre-aligned BAM or CRAM as input so a user can
run ``--stage calling`` / ``--stage filtering`` directly on an external
alignment without re-mapping. The ingest step normalizes the alignment into the
existing calling input anchor (``results/mapping/{sample}/{sample}.cram``) so
``calling.smk`` / ``filtering.smk`` need zero changes.

Test scope markers:
  * Local (default): pure-Python logic — run.py input resolution,
    make_samples_tsv bam-row builder + schema-agnostic write_tsv,
    auto_params bam-mode, and the sequence-dictionary validation helper.
  * Cluster-only (``integration`` marker): ``snakemake -n`` dry-run that
    requires the snakemake binary; skipped automatically when unavailable.
"""

from __future__ import annotations

import csv
import importlib.util
import shutil
from pathlib import Path
from typing import Any

import pytest
import yaml
from scripts import auto_params
from scripts import make_samples_tsv as mst
from scripts import validate_alignment_ref as var

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def _load_run_module() -> Any:
    """Load ``run.py`` as a module without executing ``main()``."""
    spec = importlib.util.spec_from_file_location("bsmn_run", PROJECT_ROOT / "run.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _touch(path: Path) -> None:
    """Create an empty placeholder file (content is irrelevant for path tests)."""
    path.write_bytes(b"")


# ---------------------------------------------------------------------------
# 1. run.py _resolve_inputs — BAM-mode detection
# ---------------------------------------------------------------------------


def test_resolve_inputs_single_bam_returns_bam_row(tmp_path: Path) -> None:
    """A single ``.bam`` path yields one bam-mode row keyed sample_id/bam."""
    bam = tmp_path / "HG002.bam"
    _touch(bam)

    run_mod = _load_run_module()
    rows = run_mod._resolve_inputs([str(bam)], recursive=False, pattern=None)

    assert len(rows) == 1
    row = rows[0]
    assert set(row.keys()) == {"sample_id", "bam"}
    assert row["sample_id"] == "HG002"
    assert row["bam"] == str(bam.resolve())


def test_resolve_inputs_single_cram_returns_bam_row(tmp_path: Path) -> None:
    """A single ``.cram`` path also yields a bam-mode row (bam value = cram path)."""
    cram = tmp_path / "HG002.cram"
    _touch(cram)

    run_mod = _load_run_module()
    rows = run_mod._resolve_inputs([str(cram)], recursive=False, pattern=None)

    assert len(rows) == 1
    assert rows[0]["sample_id"] == "HG002"
    assert rows[0]["bam"] == str(cram.resolve())


def test_resolve_inputs_bam_sample_id_uses_full_stem_without_pattern(tmp_path: Path) -> None:
    """Without --pattern, the full filename stem is the sample_id."""
    bam = tmp_path / "HG002.GRCh38.300x_chr20.bam"
    _touch(bam)

    run_mod = _load_run_module()
    rows = run_mod._resolve_inputs([str(bam)], recursive=False, pattern=None)

    assert rows[0]["sample_id"] == "HG002.GRCh38.300x_chr20"


def test_resolve_inputs_bam_sample_id_applies_pattern(tmp_path: Path) -> None:
    """With --pattern, the regex named group 'sample_id' is extracted from the stem."""
    bam = tmp_path / "HG002.GRCh38.300x_chr20.bam"
    _touch(bam)

    run_mod = _load_run_module()
    rows = run_mod._resolve_inputs(
        [str(bam)], recursive=False, pattern=r"(?P<sample_id>[^.]+)"
    )

    assert rows[0]["sample_id"] == "HG002"


def test_resolve_inputs_directory_of_bams(tmp_path: Path) -> None:
    """A directory containing alignment files (and no FASTQ) yields one row each."""
    _touch(tmp_path / "A.bam")
    _touch(tmp_path / "B.cram")

    run_mod = _load_run_module()
    rows = run_mod._resolve_inputs([str(tmp_path)], recursive=False, pattern=None)

    by_id = {r["sample_id"]: r for r in rows}
    assert set(by_id) == {"A", "B"}
    assert all(set(r.keys()) == {"sample_id", "bam"} for r in rows)


# ---------------------------------------------------------------------------
# 1b. run.py _resolve_inputs — FASTQ regression (must be unchanged)
# ---------------------------------------------------------------------------


def _make_gz(path: Path) -> None:
    import gzip

    with gzip.open(path, "wt") as fh:
        fh.write("@r1\nACGT\n+\nIIII\n")


def test_resolve_inputs_fastq_pair_still_returns_fastq_row(tmp_path: Path) -> None:
    """Regression: a single R1 FASTQ still produces a fastq-schema row."""
    r1 = tmp_path / "sample_R1.fastq.gz"
    r2 = tmp_path / "sample_R2.fastq.gz"
    _make_gz(r1)
    _make_gz(r2)

    run_mod = _load_run_module()
    rows = run_mod._resolve_inputs([str(r1)], recursive=False, pattern=None)

    assert len(rows) == 1
    assert set(rows[0].keys()) == {"sample_id", "readgroup", "fq1", "fq2"}
    assert rows[0]["sample_id"] == "sample"


def test_resolve_inputs_directory_of_fastqs_unchanged(tmp_path: Path) -> None:
    """Regression: a directory of FASTQ pairs still scans into fastq rows."""
    _make_gz(tmp_path / "s1_R1.fastq.gz")
    _make_gz(tmp_path / "s1_R2.fastq.gz")

    run_mod = _load_run_module()
    rows = run_mod._resolve_inputs([str(tmp_path)], recursive=False, pattern=None)

    assert len(rows) == 1
    assert "fq1" in rows[0]
    assert "bam" not in rows[0]


# ---------------------------------------------------------------------------
# 2. make_samples_tsv — bam-row builder + schema-agnostic write_tsv
# ---------------------------------------------------------------------------


def test_build_bam_row_full_stem(tmp_path: Path) -> None:
    bam = tmp_path / "HG002.GRCh38.300x_chr20.bam"
    _touch(bam)

    row = mst.build_bam_row(bam, pattern=None)
    assert row == {"sample_id": "HG002.GRCh38.300x_chr20", "bam": str(bam.resolve())}


def test_build_bam_row_with_pattern(tmp_path: Path) -> None:
    cram = tmp_path / "HG002.GRCh38.cram"
    _touch(cram)

    row = mst.build_bam_row(cram, pattern=r"(?P<sample_id>[^.]+)")
    assert row["sample_id"] == "HG002"
    assert row["bam"] == str(cram.resolve())


def test_write_tsv_bam_schema_roundtrip(tmp_path: Path) -> None:
    """write_tsv writes the bam schema and a DictReader round-trips it."""
    out = tmp_path / "samples.tsv"
    rows = [{"sample_id": "HG002", "bam": "/data/HG002.bam"}]

    mst.write_tsv(rows, str(out))

    text = out.read_text(encoding="utf-8")
    assert text.splitlines()[0] == "sample_id\tbam"

    with open(out) as fh:
        read_back = list(csv.DictReader(fh, delimiter="\t"))
    assert read_back == rows


def test_write_tsv_fastq_schema_still_works(tmp_path: Path) -> None:
    """Regression: the fastq schema header/columns are unchanged."""
    out = tmp_path / "samples.tsv"
    rows = [{"sample_id": "S1", "readgroup": "RG1", "fq1": "/a/1.fq.gz", "fq2": "/a/2.fq.gz"}]

    mst.write_tsv(rows, str(out))

    assert out.read_text(encoding="utf-8").splitlines()[0] == "sample_id\treadgroup\tfq1\tfq2"


def test_write_tsv_keeps_sample_id_first(tmp_path: Path) -> None:
    """Column order is sample_id-first regardless of dict insertion order."""
    out = tmp_path / "samples.tsv"
    rows = [{"bam": "/data/x.bam", "sample_id": "X"}]

    mst.write_tsv(rows, str(out))

    assert out.read_text(encoding="utf-8").splitlines()[0] == "sample_id\tbam"


# ---------------------------------------------------------------------------
# 3. auto_params — bam-mode (skip FASTQ opening, emit resource keys)
# ---------------------------------------------------------------------------


def test_auto_params_bam_mode_emits_required_keys(tmp_path: Path) -> None:
    """A bam-schema samples.tsv produces resolved_params with all resource keys
    consumed by calling.smk / filtering.smk / mapping.smk, and NO FASTQ access."""
    samples = tmp_path / "samples.tsv"
    samples.write_text("sample_id\tbam\nHG002\t/data/HG002.bam\n", encoding="utf-8")
    out = tmp_path / "resolved_params.yaml"

    params = auto_params.resolve_params_from_samples(str(samples), str(out))

    # Keys consumed by the rules (RESOLVED.get / RESOLVED[...]):
    for key in ("gatk_memory_gb", "bqsr_memory_gb", "markdup_memory"):
        assert key in params, f"missing required resource key: {key}"

    assert params["sequencer"] == "external-alignment"
    assert params["optical_duplicate_pixel_distance"] is None

    # The written YAML must contain the same keys.
    written = yaml.safe_load(out.read_text(encoding="utf-8"))
    assert written["sequencer"] == "external-alignment"
    assert written["gatk_memory_gb"] == params["gatk_memory_gb"]


def test_auto_params_bam_mode_does_not_open_alignment(tmp_path: Path, monkeypatch) -> None:
    """The bam path must never be opened as a FASTQ (no header read)."""
    samples = tmp_path / "samples.tsv"
    samples.write_text("sample_id\tbam\nHG002\t/data/does_not_exist.bam\n", encoding="utf-8")
    out = tmp_path / "resolved_params.yaml"

    def _boom(_path: str) -> str:
        raise AssertionError("FASTQ header was read in bam-mode")

    monkeypatch.setattr(auto_params, "_read_first_header", _boom)

    # Must succeed without touching the (nonexistent) bam as a FASTQ.
    auto_params.resolve_params_from_samples(str(samples), str(out))
    assert out.exists()


def test_auto_params_fastq_mode_unchanged(tmp_path: Path) -> None:
    """Regression: a fastq-schema samples.tsv still detects the sequencer."""
    r1 = tmp_path / "r1.fastq"
    # NovaSeq 6000 instrument id (A00xxx) → patterned flowcell.
    r1.write_text("@A00100:1:HXXX:1:1101:1000:1000 1:N:0:1\nACGT\n+\nIIII\n", encoding="utf-8")
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        f"sample_id\treadgroup\tfq1\tfq2\nS1\tRG1\t{r1}\t{r1}\n", encoding="utf-8"
    )
    out = tmp_path / "resolved_params.yaml"

    params = auto_params.resolve_params_from_samples(str(samples), str(out))
    assert params["sequencer"] == "NovaSeq 6000"
    assert params["optical_duplicate_pixel_distance"] == 2500


# ---------------------------------------------------------------------------
# 4. validate_alignment_ref — sequence-dictionary compatibility logic
# ---------------------------------------------------------------------------


def _write_dict(path: Path, contigs: list[tuple[str, int]]) -> None:
    """Write a minimal Picard .dict (@HD + @SQ lines)."""
    lines = ["@HD\tVN:1.6"]
    for name, ln in contigs:
        lines.append(f"@SQ\tSN:{name}\tLN:{ln}\tUR:file:/ref.fasta")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_parse_dict_contigs(tmp_path: Path) -> None:
    d = tmp_path / "ref.dict"
    _write_dict(d, [("chr1", 248956422), ("chr20", 64444167)])
    contigs = var.parse_dict_contigs(d)
    assert contigs == {"chr1": 248956422, "chr20": 64444167}


def test_parse_sq_lines() -> None:
    header = "\n".join(
        [
            "@HD\tVN:1.6\tSO:coordinate",
            "@SQ\tSN:chr1\tLN:248956422",
            "@SQ\tSN:chr20\tLN:64444167",
            "@RG\tID:x\tSM:HG002",
        ]
    )
    contigs = var.parse_sq_lines(header)
    assert contigs == {"chr1": 248956422, "chr20": 64444167}


def test_compatibility_subset_is_ok() -> None:
    """A BAM with a subset of reference contigs (matching lengths) is compatible."""
    ref = {"chr1": 248956422, "chr20": 64444167, "chrX": 156040895}
    aln = {"chr20": 64444167}  # subset
    errors = var.check_compatibility(aln, ref)
    assert errors == []


def test_compatibility_extra_decoy_contigs_ok() -> None:
    """Extra decoy/alt contigs only in the alignment do not fail (warn, not error)."""
    ref = {"chr1": 248956422}
    aln = {"chr1": 248956422, "chrUn_decoy1": 1000}
    errors = var.check_compatibility(aln, ref)
    assert errors == []


def test_compatibility_length_mismatch_fails() -> None:
    """A shared contig name with a different length is a hard failure naming the contig."""
    ref = {"chr1": 248956422}
    aln = {"chr1": 249250621}  # GRCh37 chr1 length — incompatible
    errors = var.check_compatibility(aln, ref)
    assert len(errors) == 1
    assert "chr1" in errors[0]
    assert "248956422" in errors[0]
    assert "249250621" in errors[0]


def test_compatibility_no_shared_contigs_fails() -> None:
    """Zero overlapping contig names means the references are unrelated → fail."""
    ref = {"chr1": 248956422}
    aln = {"1": 248956422}  # no 'chr' prefix → no name overlap
    errors = var.check_compatibility(aln, ref)
    assert errors
    assert any("no shared contigs" in e.lower() for e in errors)


def test_validate_against_dict_passes(tmp_path: Path) -> None:
    d = tmp_path / "ref.dict"
    _write_dict(d, [("chr1", 248956422), ("chr20", 64444167)])
    header = "@HD\tVN:1.6\n@SQ\tSN:chr20\tLN:64444167\n"
    errors = var.validate_against_dict(header, d)
    assert errors == []


def test_validate_against_dict_reports_mismatch(tmp_path: Path) -> None:
    d = tmp_path / "ref.dict"
    _write_dict(d, [("chr20", 64444167)])
    header = "@HD\tVN:1.6\n@SQ\tSN:chr20\tLN:99999999\n"
    errors = var.validate_against_dict(header, d)
    assert errors
    assert "chr20" in errors[0]


# ---------------------------------------------------------------------------
# 5. Cluster-only: snakemake -n dry-run on a tiny fixture
# ---------------------------------------------------------------------------


@pytest.mark.integration
@pytest.mark.skipif(
    shutil.which("snakemake") is None,
    reason="cluster-only: requires the snakemake binary (not installed locally)",
)
def test_ingest_dryrun_selects_ingest_not_mapping(tmp_path: Path) -> None:
    """A bam-mode samples.tsv must produce a DAG containing ingest_alignment +
    mutect2_scatter and NO mapping rules. Runs only where snakemake exists."""
    import os
    import subprocess

    workdir = tmp_path
    (workdir / "config").mkdir()
    (workdir / "resources" / "hg38").mkdir(parents=True)

    # Minimal fake reference + dict (chr20 only, matching the fixture bam header).
    ref_fasta = workdir / "resources" / "hg38" / "ref.fasta"
    ref_fasta.write_text(">chr20\nACGT\n", encoding="utf-8")
    ref_dict = workdir / "resources" / "hg38" / "ref.dict"
    _write_dict(ref_dict, [("chr20", 64444167)])

    fake_bam = workdir / "HG002.bam"
    _touch(fake_bam)

    # samples.tsv (bam schema) + resolved_params.yaml (skip auto_params subprocess).
    (workdir / "config" / "samples.tsv").write_text(
        f"sample_id\tbam\nHG002\t{fake_bam}\n", encoding="utf-8"
    )
    var_out = workdir / "config" / "resolved_params.yaml"
    auto_params.resolve_params_from_samples(
        str(workdir / "config" / "samples.tsv"), str(var_out)
    )

    # containers.yaml — only the keys the included rules reference.
    shutil.copy(PROJECT_ROOT / "config" / "containers.yaml", workdir / "config" / "containers.yaml")

    # chr20 config pointing at the fake reference.
    cfg = yaml.safe_load((PROJECT_ROOT / "config" / "config.chr20.yaml").read_text())
    cfg["ref"]["fasta"] = "resources/hg38/ref.fasta"
    cfg["ref"]["dict"] = "resources/hg38/ref.dict"
    (workdir / "config" / "config.chr20.yaml").write_text(yaml.dump(cfg), encoding="utf-8")

    proc = subprocess.run(
        [
            "snakemake",
            "--snakefile",
            str(PROJECT_ROOT / "workflow" / "Snakefile"),
            "--configfile",
            "config/config.chr20.yaml",
            "-n",
            "-p",
        ],
        cwd=workdir,
        capture_output=True,
        text=True,
        env={**os.environ},
    )
    out = proc.stdout + proc.stderr
    assert "ingest_alignment" in out, out
    assert "mutect2_scatter" in out, out
    assert "bwa_mem" not in out and "apply_bqsr" not in out, out
