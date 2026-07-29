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
# Host-aware thread / engine derivation  (the auto-tuner)
#
# @MX:NOTE: [AUTO] These pick each job's WIDTH (thread count) and the
#           MarkDuplicates engine from the host's core budget + RAM, so a run needs
#           no hand-written per-host overlay (e.g. config.storage32.yaml). How MANY
#           jobs run at once is still governed by the memory admission budget above
#           (concurrency = --resources mem_mb // per-job mem_mb). Evidence tier per
#           knob, verified against primary sources:
#             [DOC]       markdup_threads ceiling 16 — MarkDuplicatesSpark scales
#                         linearly only to ~16 cores (GATK Javadoc); mutect2 threads
#                         default 4 and PairHMM-only (GATK source); BQSR takes no
#                         thread knob (Broad WDL — parallelised by interval scatter).
#             [HEURISTIC] bwa per-job width and budgeting on the logical core count
#                         — no published optimum exists, so these are labelled and
#                         overridable via config ('auto' defers here; an int/str
#                         pins a value). markdup defaults to Picard; spark is an
#                         explicit opt-in (needs fast local scratch — see FINDING in
#                         the auto-tuner memory / derive_markdup docstring).
# @MX:REASON: fan_in >= 3 — emitted into resolved_params.yaml and read by
#             mapping.smk (bwa/markdup) and calling.smk (mutect2).
# ---------------------------------------------------------------------------

# [DOC] MarkDuplicatesSpark scales linearly only to ~16 cores — never exceed.
_MARKDUP_SPARK_MAX_THREADS = 16
# [DOC] Mutect2 --native-pair-hmm-threads default; only the PairHMM is threaded, so
# surplus cores are better spent on concurrent chromosomes than on this knob.
_MUTECT2_PAIRHMM_THREADS = 4
# [HEURISTIC] Preferred per-job bwa-mem width when cores are plentiful. bwa-mem has
# no published thread optimum, so this only *targets* a moderate width and lets
# extra cores buy more concurrent read-groups. Benchmark per host-class to confirm.
_BWA_TARGET_THREADS = 12
# One bwa_mem_sort job's working set: hg38 BWT index (~5.5 GB) + sort buffer
# (~8 GB) + overhead. Bounds bwa concurrency by RAM so cores are not
# oversubscribed against memory on RAM-poor / many-core hosts.
_BWA_JOB_GB = 15


def derive_bwa_threads(cores: int, ram_gb: int) -> int:
    """Per-job bwa-mem thread count  ([HEURISTIC] — no published optimum).

    Concurrency-aware: aim for ~``_BWA_TARGET_THREADS``-wide jobs, cap the number of
    concurrent jobs by what RAM allows (~``_BWA_JOB_GB`` each), then give each job an
    equal share of the cores. Small hosts collapse to one job using every core;
    RAM-bound many-core hosts get fewer, wider jobs so cores are not left idle behind
    a memory ceiling. Override via ``mapping.bwa_threads`` in config.
    """
    cores = max(1, int(cores))
    ram_gb = max(1, int(ram_gb))
    ram_concurrency = max(1, ram_gb // _BWA_JOB_GB)
    core_concurrency = max(1, cores // _BWA_TARGET_THREADS)
    concurrency = max(1, min(core_concurrency, ram_concurrency))
    return max(1, cores // concurrency)


def derive_markdup(cores: int) -> tuple[str, int]:
    """(engine, threads) for MarkDuplicates.

    Auto-default is **Picard**, unconditionally. On a multi-sample cohort Picard
    pipelines across samples on a many-core host (measured: 22/30 markdups flowed
    through on the 168-core VM without stalling), needs no fast scratch, and has no
    ``.parts`` staging fragility. Spark is left as an explicit opt-in
    (``markdup_engine: spark`` in config) because it only wins with fast LOCAL
    scratch — measured ~2x SLOWER than Picard when its shuffle spills to NFS — and
    warrants per-host validation. The returned thread count (capped at 16, where
    MarkDuplicatesSpark scaling flattens per [DOC]) applies only when spark is
    pinned; Picard ignores it.
    """
    cores = max(1, int(cores))
    return "picard", min(_MARKDUP_SPARK_MAX_THREADS, cores)


def derive_mutect2_threads(cores: int) -> int:
    """--native-pair-hmm-threads  ([DOC] GATK default 4).

    Only the PairHMM is threaded, so raising this barely moves wall-clock; surplus
    cores go to concurrent chromosomes instead.
    """
    return min(_MUTECT2_PAIRHMM_THREADS, max(1, int(cores)))


def get_logical_cpus() -> int:
    """Logical CPU count (hyperthreads included) — the budget ``--cores`` schedules
    against."""
    return max(1, os.cpu_count() or 1)


def get_physical_cpus() -> int:
    """Physical core count; falls back to the logical count when undeterminable.
    Reported for transparency; the tuner budgets on the logical count ([HEURISTIC]
    — no bwa/GATK authority prescribes physical-vs-logical)."""
    try:
        n = psutil.cpu_count(logical=False)
        return int(n) if n else get_logical_cpus()
    except Exception:
        return get_logical_cpus()


def get_total_memory_gb() -> int:
    """Return total system memory in GB (integer, floored)."""
    try:
        return int(psutil.virtual_memory().total >> 30)
    except Exception:
        return 16


def _tuning_params(cores: int | None = None) -> dict[str, Any]:
    """Host-derived resource + thread/engine knobs consumed by the workflow rules.

    ``cores`` is the budget to tune per-job widths for (default: detected logical
    CPUs) — pass the run's ``--cores`` to match the real budget. Keys read by rules:
      * mapping.smk   : bwa_threads, sort_threads, sort_memory, markdup_engine,
                        markdup_threads, markdup_memory, bqsr_memory_gb
      * calling.smk   : mutect2_threads, gatk_memory_gb, bqsr_memory_gb
      * filtering.smk : gatk_memory_gb
    """
    n_cores = int(cores) if cores else get_logical_cpus()
    ram_gb = get_total_memory_gb()
    bwa = derive_bwa_threads(n_cores, ram_gb)
    markdup_engine, markdup_threads = derive_markdup(n_cores)
    return {
        "system_total_memory_mb": get_total_memory_mb(),
        "system_cpus": get_logical_cpus(),
        "system_physical_cpus": get_physical_cpus(),
        "tuned_for_cores": n_cores,
        "bwa_threads": bwa,
        "sort_threads": min(4, bwa),
        "sort_memory": get_sort_memory(),
        "markdup_engine": markdup_engine,
        "markdup_threads": markdup_threads,
        "mutect2_threads": derive_mutect2_threads(n_cores),
        "bqsr_memory_gb": get_bqsr_memory_gb(),
        "markdup_memory": f"{get_markdup_memory_gb()}G",
        "gatk_memory_gb": get_gatk_memory_gb(),
    }


# ---------------------------------------------------------------------------
# Main resolver
# ---------------------------------------------------------------------------


def resolve_params(
    fastq_path: str, output_yaml: str, cores: int | None = None
) -> dict[str, Any]:
    """
    Analyze FASTQ, determine all pipeline parameters, and write
    resolved_params.yaml.

    ``cores`` tunes the per-job thread widths (default: detected logical CPUs).
    Returns the parameter dict (includes sequencer_evidence).
    """
    sequencer, evidence = _detect_sequencer_with_evidence(fastq_path)

    params: dict[str, Any] = {
        "resolved_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "input_fastq": os.path.abspath(fastq_path),
        "sequencer": sequencer,
        "sequencer_evidence": evidence,
        "optical_duplicate_pixel_distance": get_optical_duplicate_pixel_distance(sequencer),
        **_tuning_params(cores),
    }

    out_dir = os.path.dirname(output_yaml)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(output_yaml, "w") as fh:
        yaml.dump(params, fh, default_flow_style=False, sort_keys=False)

    return params


def resolve_params_bam(
    alignment_path: str, output_yaml: str, cores: int | None = None
) -> dict[str, Any]:
    """Resolve parameters for an external-alignment (BAM/CRAM) input.

    Skips FASTQ opening and sequencer detection entirely: a pre-aligned input
    has no FASTQ header to read. Emits the same tuning keys the rules consume
    so calling/filtering run unchanged, with ``sequencer`` set to a sentinel and
    ``optical_duplicate_pixel_distance`` null (no de-duplication is performed in
    ingest mode). ``cores`` tunes per-job thread widths.
    """
    params: dict[str, Any] = {
        "resolved_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "input_alignment": os.path.abspath(alignment_path),
        "schema": "bam",
        "sequencer": "external-alignment",
        "optical_duplicate_pixel_distance": None,
        **_tuning_params(cores),
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
def resolve_params_from_samples(
    samples_tsv: str, output_yaml: str, cores: int | None = None
) -> dict[str, Any]:
    """Detect the samples.tsv schema and resolve parameters accordingly.

    * bam schema   → :func:`resolve_params_bam` (no FASTQ access)
    * fastq schema → :func:`resolve_params` on the first ``fq1`` entry

    ``cores`` is forwarded to tune per-job thread widths for the run's budget.
    """
    with open(samples_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        first = next(reader)

    if _detect_schema(header) == "bam":
        return resolve_params_bam(first["bam"], output_yaml, cores=cores)
    return resolve_params(first["fq1"], output_yaml, cores=cores)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def _print_tuning(params: dict[str, Any]) -> None:
    """Print the host-tuned thread/engine choices with their evidence tier, so the
    user can see (and, if needed, override) exactly what the auto-tuner picked."""
    ram_gb = params.get("system_total_memory_mb", 0) // 1024
    print(
        f"[auto_params] host             : {params['system_cpus']} logical / "
        f"{params['system_physical_cpus']} physical cores, {ram_gb} GB RAM "
        f"(tuned for {params['tuned_for_cores']} cores)"
    )
    print(
        f"[auto_params] bwa_threads      : {params['bwa_threads']}"
        "        [HEURISTIC — no published optimum; benchmark to confirm]"
    )
    print(
        f"[auto_params] markdup          : {params['markdup_engine']} x "
        f"{params['markdup_threads']}   [DOC: Spark linear to ~16 cores]"
    )
    print(
        f"[auto_params] mutect2_threads  : {params['mutect2_threads']}"
        "        [DOC: GATK default 4, PairHMM-only]"
    )


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
    parser.add_argument(
        "--cores",
        type=int,
        default=None,
        help="Core budget to tune per-job thread widths for "
        "(default: detected logical CPUs). Pass the run's --cores to match it.",
    )
    args = parser.parse_args()

    if args.samples:
        params = resolve_params_from_samples(args.samples, args.output, cores=args.cores)
    else:
        params = resolve_params(args.fastq, args.output, cores=args.cores)

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
        _print_tuning(params)
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
    _print_tuning(params)
    print(f"[auto_params] saved to         → {args.output}")


if __name__ == "__main__":
    main()
