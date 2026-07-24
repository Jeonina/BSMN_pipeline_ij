#!/usr/bin/env python3
"""
auto_params.py — Auto-determine BSMN pipeline parameters from FASTQ input.

Detects sequencer type from instrument ID pattern, colon-field count, and
flowcell ID; queries system resources; writes config/resolved_params.yaml
for reproducible runs.

Usage:
    python scripts/auto_params.py <R1_fastq>
    python scripts/auto_params.py --samples config/samples.tsv
    python scripts/auto_params.py <R1_fastq> --output config/resolved_params.yaml

Sequencer detection priority (PIPELINE.md):
    1. Instrument ID regex  — most reliable
    2. Colon-field count    — >= 8 fields → patterned flowcell heuristic
    3. Unknown fallback     — ODPD=100, conservative

OPTICAL_DUPLICATE_PIXEL_DISTANCE (ODPD):
    Patterned flowcell   : NovaSeq 6000, NovaSeq X, HiSeq X          → 2500
    Unpatterned flowcell : HiSeq 2500/3000/4000, MiSeq, NextSeq, Unknown → 100
"""

import argparse
import csv
import datetime
import gzip
import os
import re
from collections.abc import Sequence
from typing import Any

import psutil
import yaml

# ---------------------------------------------------------------------------
# Sequencer identification tables
# ---------------------------------------------------------------------------

# Patterned flowcell sequencers → ODPD = 2500
# Note: HiSeq 3000/4000 use patterned hardware but BSMN pipeline treats them
# as unpatterned per lab convention (ODPD=100) to avoid over-deduplication.
_PATTERNED_SEQUENCERS: frozenset[str] = frozenset(
    {
        "NovaSeq 6000",
        "NovaSeq X",
        "HiSeq X",
        "NovaSeq",  # legacy name — kept for backward compatibility
    }
)

# Instrument ID prefix → sequencer name
# Ordered most-specific first to avoid false matches.
# Sources: Illumina BaseSpace naming, ENCODE portal metadata, SRA instrument fields
_INSTRUMENT_RULES: list[tuple[re.Pattern[str], str]] = [
    (re.compile(r"^LH\d"), "NovaSeq X"),  # e.g. LH00204, LH00478
    (re.compile(r"^A\d{5}$"), "NovaSeq 6000"),  # e.g. A00100, A00266
    (re.compile(r"^[EK]\d{5}$"), "HiSeq X"),  # e.g. E00143, K00145
    (re.compile(r"^J\d{5}$"), "HiSeq 3000/4000"),  # e.g. J00120, J00122
    (re.compile(r"^(SN|D)\d+"), "HiSeq 2500"),  # e.g. SN0196, D00195
    (re.compile(r"^HWI"), "HiSeq 2500"),  # e.g. HWI-ST*, HWI-M*
    (re.compile(r"^M\d{5}$"), "MiSeq"),  # e.g. M03213
    (re.compile(r"^(NS|NB)\d+"), "NextSeq 500/550"),  # e.g. NS500487, NB501234
    (re.compile(r"^VH\d+"), "NextSeq 2000"),  # e.g. VH00204
]


# ---------------------------------------------------------------------------
# Sequencer detection (internal)
# ---------------------------------------------------------------------------


def _read_first_header(fastq_path: str) -> str:
    """Read the FASTQ header line from a gzipped or plain file."""
    opener = gzip.open if fastq_path.endswith(".gz") else open
    with opener(fastq_path, "rt") as fh:
        return fh.readline().strip()


def _detect_sequencer_with_evidence(
    fastq_path: str,
) -> tuple[str, dict[str, Any]]:
    """
    Return (sequencer_name, evidence_dict).

    Evidence keys: instrument_id, colon_fields, flowcell_id,
                   detection_method, matched_pattern/note, read_name_example.
    """
    header = _read_first_header(fastq_path)
    read_name = header.lstrip("@").split()[0]
    fields = read_name.split(":")
    n_colons = read_name.count(":")
    instrument = fields[0]
    # CASAVA 1.8+ layout: instrument:run:flowcell:lane:tile:x:y
    flowcell = fields[2] if len(fields) >= 3 else ""

    base: dict[str, Any] = {
        "instrument_id": instrument,
        "colon_fields": n_colons + 1,
        "flowcell_id": flowcell,
        "read_name_example": read_name,
    }

    # Priority 1 — instrument ID pattern
    for pattern, name in _INSTRUMENT_RULES:
        if pattern.match(instrument):
            return name, {
                **base,
                "detection_method": "instrument_id_pattern",
                "matched_pattern": pattern.pattern,
            }

    # Priority 2 — colon count heuristic
    # CASAVA 1.8+ base format = 7 fields (6 colons).
    # Patterned-flowcell runs commonly add a UMI field → 8 fields (7 colons).
    if n_colons >= 7:
        return "NovaSeq 6000", {
            **base,
            "detection_method": "colon_count_heuristic",
            "note": (
                f"Instrument '{instrument}' not in known ID list; "
                f"{n_colons + 1} colon-fields suggest patterned flowcell (NovaSeq assumed)"
            ),
        }

    # Priority 3 — unknown fallback (conservative ODPD=100)
    return "Unknown", {
        **base,
        "detection_method": "unknown_fallback",
        "note": (
            f"Instrument '{instrument}' not recognized; defaulting to ODPD=100 (conservative)"
        ),
    }


# ---------------------------------------------------------------------------
# Public sequencer API
# ---------------------------------------------------------------------------


def detect_sequencer(fastq_path: str) -> str:
    """
    Return sequencer name from the first FASTQ read name.

    Detection priority:
      1. Instrument ID pattern (e.g. A00100 → NovaSeq 6000)
      2. Colon-field count   (>= 8 fields → NovaSeq 6000 assumed)
      3. Unknown fallback
    """
    name, _ = _detect_sequencer_with_evidence(fastq_path)
    return name


def get_optical_duplicate_pixel_distance(sequencer: str) -> int:
    """
    Return Picard OPTICAL_DUPLICATE_PIXEL_DISTANCE for the sequencer.

    Patterned flowcell (NovaSeq 6000, NovaSeq X, HiSeq X): 2500
    All others (HiSeq 2500/3000/4000, MiSeq, NextSeq, Unknown): 100
    """
    return 2500 if sequencer in _PATTERNED_SEQUENCERS else 100


# ---------------------------------------------------------------------------
# System resource queries
# ---------------------------------------------------------------------------


def get_bwa_threads() -> int:
    """Return CPU count capped at 4 on low-memory systems (< 8 GB)."""
    cpus = max(1, os.cpu_count() or 1)
    try:
        total_gb: int = psutil.virtual_memory().total >> 30
        if total_gb < 8:
            return min(cpus, 4)
    except Exception:
        pass
    return int(cpus)


# ---------------------------------------------------------------------------
# Per-job memory footprints
#
# @MX:NOTE: [AUTO] These are per-JOB working-set sizes, deliberately NOT scaled
#           to total system RAM. Each value maps 1:1 onto a rule's
#           `resources.mem_mb`, which Snakemake uses as an admission-control
#           budget: concurrency = (--resources mem_mb) // (per-job mem_mb).
#           Deriving them from total RAM therefore INVERTS the intended effect —
#           a bigger machine reserved a bigger heap per job and ran FEWER jobs.
#           Measured on a 180-core/98 GB host: total-RAM scaling produced 24-49 GB
#           reservations against 13 GB actual RSS, pinning the per-chromosome
#           scatter rules (base_recalibrator, mutect2) to 1-3 concurrent jobs and
#           leaving 177 of 180 cores idle. More RAM must buy more CONCURRENCY,
#           not a larger heap. Values below are the tools' real requirements.
# @MX:REASON: fan_in >= 3 — read by mapping.smk, calling.smk, and filtering.smk
#             via config/resolved_params.yaml.
# ---------------------------------------------------------------------------

# GATK streaming walkers (BaseRecalibrator, ApplyBQSR, Mutect2, pileups).
# These stream over the alignment and hold a bounded window in memory.
_GATK_HEAP_GB = 8

# Picard/Spark MarkDuplicates keeps a read-name → position map for unpaired and
# not-yet-mated reads, so it needs materially more than the streaming walkers.
_MARKDUP_HEAP_GB = 16

# sambamba sort in-memory buffer before spilling to disk (per bwa_mem_sort job).
_SORT_MEMORY_GB = 8

# Never let a single job's heap exceed this fraction of system RAM, so the
# defaults degrade gracefully on small hosts (e.g. a 32 GB workstation).
_MAX_HEAP_FRACTION = 4


def _cap_to_system(heap_gb: int, floor_gb: int = 2) -> int:
    """Clamp a per-job heap to ``1/_MAX_HEAP_FRACTION`` of total system RAM.

    Returns ``heap_gb`` unchanged on hosts with ample RAM; shrinks it (never
    below ``floor_gb``) on small hosts so a job is still schedulable.
    """
    try:
        total_gb: int = psutil.virtual_memory().total >> 30
    except Exception:
        return heap_gb
    return max(floor_gb, min(heap_gb, total_gb // _MAX_HEAP_FRACTION))


def get_sort_memory() -> str:
    """Return the sambamba sort buffer for one bwa_mem_sort job.

    Fixed rather than a share of *available* RAM: the previous
    30%-of-available rule made the resolved parameters depend on whatever else
    happened to be running at resolve time, so two runs of the same cohort
    could sort with different buffer sizes.
    """
    return f"{_cap_to_system(_SORT_MEMORY_GB)}GB"


def get_bqsr_memory_gb() -> int:
    """Return the Java heap for GATK streaming walkers (BQSR, Mutect2)."""
    return _cap_to_system(_GATK_HEAP_GB, floor_gb=4)


def get_markdup_memory_gb() -> int:
    """Return the Java heap for MarkDuplicates (Picard or Spark)."""
    return _cap_to_system(_MARKDUP_HEAP_GB, floor_gb=4)


def get_gatk_memory_gb() -> int:
    """Return the Java heap for general GATK tools (MergeVcfs, gathers, etc.)."""
    return _cap_to_system(_GATK_HEAP_GB, floor_gb=4)


def get_total_memory_mb() -> int:
    """Return total system memory in MB."""
    try:
        return int(psutil.virtual_memory().total >> 20)
    except Exception:
        return 16384


# ---------------------------------------------------------------------------
# Main resolver
# ---------------------------------------------------------------------------


def resolve_params(fastq_path: str, output_yaml: str) -> dict[str, Any]:
    """
    Analyze FASTQ, determine all pipeline parameters, and write
    resolved_params.yaml.

    Returns the parameter dict (includes sequencer_evidence).
    """
    sequencer, evidence = _detect_sequencer_with_evidence(fastq_path)
    bqsr_mem = get_bqsr_memory_gb()
    markdup_mem = get_markdup_memory_gb()

    gatk_mem = get_gatk_memory_gb()

    params: dict[str, Any] = {
        "resolved_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "input_fastq": os.path.abspath(fastq_path),
        "sequencer": sequencer,
        "sequencer_evidence": evidence,
        "system_total_memory_mb": get_total_memory_mb(),
        "system_cpus": max(1, os.cpu_count() or 1),
        "bwa_threads": get_bwa_threads(),
        "sort_threads": min(4, get_bwa_threads()),
        "sort_memory": get_sort_memory(),
        "optical_duplicate_pixel_distance": get_optical_duplicate_pixel_distance(sequencer),
        "bqsr_memory_gb": bqsr_mem,
        "markdup_memory": f"{markdup_mem}G",
        "gatk_memory_gb": gatk_mem,
    }

    out_dir = os.path.dirname(output_yaml)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(output_yaml, "w") as fh:
        yaml.dump(params, fh, default_flow_style=False, sort_keys=False)

    return params


def _resource_params() -> dict[str, Any]:
    """Return the system-derived resource keys consumed by the workflow rules.

    These are the keys the rules read from ``RESOLVED``:
      * mapping.smk : ``bqsr_memory_gb`` (also ``bwa_threads``, ``sort_*``,
        ``markdup_memory``) — only included in fastq-mode, but emitted always
        so the file is schema-complete.
      * calling.smk : ``gatk_memory_gb``, ``bqsr_memory_gb``
      * filtering.smk : ``gatk_memory_gb``
    """
    return {
        "system_total_memory_mb": get_total_memory_mb(),
        "system_cpus": max(1, os.cpu_count() or 1),
        "bwa_threads": get_bwa_threads(),
        "sort_threads": min(4, get_bwa_threads()),
        "sort_memory": get_sort_memory(),
        "bqsr_memory_gb": get_bqsr_memory_gb(),
        "markdup_memory": f"{get_markdup_memory_gb()}G",
        "gatk_memory_gb": get_gatk_memory_gb(),
    }


def resolve_params_bam(alignment_path: str, output_yaml: str) -> dict[str, Any]:
    """Resolve parameters for an external-alignment (BAM/CRAM) input.

    Skips FASTQ opening and sequencer detection entirely: a pre-aligned input
    has no FASTQ header to read. Emits the same resource keys the rules consume
    so calling/filtering run unchanged, with ``sequencer`` set to a sentinel and
    ``optical_duplicate_pixel_distance`` null (no de-duplication is performed in
    ingest mode).
    """
    params: dict[str, Any] = {
        "resolved_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "input_alignment": os.path.abspath(alignment_path),
        "schema": "bam",
        "sequencer": "external-alignment",
        "optical_duplicate_pixel_distance": None,
        **_resource_params(),
    }

    out_dir = os.path.dirname(output_yaml)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(output_yaml, "w") as fh:
        yaml.dump(params, fh, default_flow_style=False, sort_keys=False)

    return params


def _detect_schema(header: Sequence[str]) -> str:
    """Return 'bam' if the samples.tsv columns are the bam schema, else 'fastq'.

    Bam schema is identified by a ``bam`` column with no ``fq1`` column.
    """
    return "bam" if ("bam" in header and "fq1" not in header) else "fastq"


# @MX:ANCHOR: [AUTO] resolve_params_from_samples — single entry point used by
#             run.py, workflow/Snakefile (auto-gen subprocess), and tests.
# @MX:REASON: this function picks the schema-correct resolver; both the Snakefile
#             auto-gen step and run.py depend on it producing a resolved_params.yaml
#             with every key calling/filtering/mapping read. fan_in >= 3.
def resolve_params_from_samples(samples_tsv: str, output_yaml: str) -> dict[str, Any]:
    """Detect the samples.tsv schema and resolve parameters accordingly.

    * bam schema   → :func:`resolve_params_bam` (no FASTQ access)
    * fastq schema → :func:`resolve_params` on the first ``fq1`` entry
    """
    with open(samples_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        first = next(reader)

    if _detect_schema(header) == "bam":
        return resolve_params_bam(first["bam"], output_yaml)
    return resolve_params(first["fq1"], output_yaml)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Auto-determine BSMN pipeline parameters from FASTQ.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("fastq", nargs="?", help="Path to R1 FASTQ (gz or plain)")
    src.add_argument(
        "--samples",
        metavar="TSV",
        help="samples.tsv; uses first fq1 entry for sequencer detection",
    )
    parser.add_argument(
        "--output",
        default="config/resolved_params.yaml",
        help="Output YAML path  (default: config/resolved_params.yaml)",
    )
    args = parser.parse_args()

    if args.samples:
        params = resolve_params_from_samples(args.samples, args.output)
    else:
        params = resolve_params(args.fastq, args.output)

    # bam-mode (external alignment) has no sequencer evidence to report.
    if params.get("schema") == "bam":
        print("[auto_params] schema           : bam (external alignment)")
        print(f"[auto_params] input_alignment  : {params['input_alignment']}")
        print(f"[auto_params] sequencer        : {params['sequencer']}")
        print(f"[auto_params] system_cpus      : {params['system_cpus']}")
        print(f"[auto_params] system_memory_mb : {params['system_total_memory_mb']}")
        print(f"[auto_params] bqsr_memory_gb   : {params['bqsr_memory_gb']}")
        print(f"[auto_params] gatk_memory_gb   : {params['gatk_memory_gb']}")
        print(f"[auto_params] markdup_memory   : {params['markdup_memory']}")
        print(f"[auto_params] saved to         → {args.output}")
        return

    ev = params["sequencer_evidence"]

    print(f"[auto_params] sequencer        : {params['sequencer']}")
    print(f"[auto_params] detection        : {ev['detection_method']}")
    print(f"[auto_params] instrument_id    : {ev['instrument_id']}")
    print(f"[auto_params] system_cpus      : {params['system_cpus']}")
    print(f"[auto_params] system_memory_mb : {params['system_total_memory_mb']}")
    print(f"[auto_params] bwa_threads      : {params['bwa_threads']}")
    print(f"[auto_params] ODPD             : {params['optical_duplicate_pixel_distance']}")
    print(f"[auto_params] bqsr_memory_gb   : {params['bqsr_memory_gb']}")
    print(f"[auto_params] gatk_memory_gb   : {params['gatk_memory_gb']}")
    print(f"[auto_params] markdup_memory   : {params['markdup_memory']}")
    print(f"[auto_params] saved to         → {args.output}")


if __name__ == "__main__":
    main()
