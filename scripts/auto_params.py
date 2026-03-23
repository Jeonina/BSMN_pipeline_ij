#!/usr/bin/env python3
"""
auto_params.py — Auto-determine BSMN pipeline parameters from FASTQ input.

Detects sequencer type from read name colon-field count, queries system
resources, and writes config/resolved_params.yaml for reproducible runs.

Usage:
    python scripts/auto_params.py <R1_fastq>
    python scripts/auto_params.py --samples config/samples.tsv
    python scripts/auto_params.py <R1_fastq> --output config/resolved_params.yaml

Sequencer detection rule (per BSMN PIPELINE.md):
    read_name.count(':') >= 7  → NovaSeq → ODPD = 2500
    otherwise                  → HiSeq   → ODPD = 100
"""

import argparse
import csv
import datetime
import gzip
import os

import psutil
import yaml


# ---------------------------------------------------------------------------
# Sequencer detection
# ---------------------------------------------------------------------------

def detect_sequencer(fastq_path: str) -> str:
    """
    Detect sequencer type from the first FASTQ read name.

    Returns 'NovaSeq' when the read name (before whitespace) contains >= 7
    colon characters; 'HiSeq' otherwise.
    """
    opener = gzip.open if fastq_path.endswith(".gz") else open
    with opener(fastq_path, "rt") as fh:
        header = fh.readline().strip()
    read_name = header.lstrip("@").split()[0]
    return "NovaSeq" if read_name.count(":") >= 7 else "HiSeq"


def get_optical_duplicate_pixel_distance(sequencer: str) -> int:
    """
    Return Picard OPTICAL_DUPLICATE_PIXEL_DISTANCE for the given sequencer.

    NovaSeq uses patterned flow cells (2500); all others use unpatterned (100).
    """
    return 2500 if sequencer == "NovaSeq" else 100


# ---------------------------------------------------------------------------
# System resource queries
# ---------------------------------------------------------------------------

def get_bwa_threads() -> int:
    """Return available logical CPU count (minimum 1)."""
    return max(1, os.cpu_count() or 1)


def get_bqsr_memory_gb() -> int:
    """
    Return recommended GATK BQSR Java heap in GB.
    Capped at half total RAM, between 4 GB and 64 GB.
    """
    try:
        total_gb: int = psutil.virtual_memory().total >> 30
        return max(4, min(total_gb // 2, 64))
    except Exception:
        return 16


def get_markdup_memory_gb() -> int:
    """Return recommended Picard MarkDuplicates Java heap in GB."""
    try:
        total_gb: int = psutil.virtual_memory().total >> 30
        return max(4, min(total_gb // 4, 32))
    except Exception:
        return 8


# ---------------------------------------------------------------------------
# Main resolver
# ---------------------------------------------------------------------------

def resolve_params(fastq_path: str, output_yaml: str) -> dict:
    """
    Analyze FASTQ, determine all pipeline parameters, and write
    resolved_params.yaml.

    Returns the parameter dict.
    """
    sequencer = detect_sequencer(fastq_path)
    bqsr_mem = get_bqsr_memory_gb()
    markdup_mem = get_markdup_memory_gb()

    params = {
        "resolved_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "input_fastq": os.path.abspath(fastq_path),
        "sequencer": sequencer,
        "bwa_threads": get_bwa_threads(),
        "sort_threads": 4,
        "sort_memory": "6GB",
        "optical_duplicate_pixel_distance": get_optical_duplicate_pixel_distance(sequencer),
        "bqsr_memory_gb": bqsr_mem,
        "markdup_memory": f"{markdup_mem}G",
    }

    out_dir = os.path.dirname(output_yaml)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(output_yaml, "w") as fh:
        yaml.dump(params, fh, default_flow_style=False, sort_keys=False)

    return params


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
        with open(args.samples) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            first = next(reader)
        fastq_path = first["fq1"]
    else:
        fastq_path = args.fastq

    params = resolve_params(fastq_path, args.output)

    print(f"[auto_params] sequencer      : {params['sequencer']}")
    print(f"[auto_params] bwa_threads    : {params['bwa_threads']}")
    print(f"[auto_params] ODPD           : {params['optical_duplicate_pixel_distance']}")
    print(f"[auto_params] bqsr_memory_gb : {params['bqsr_memory_gb']}")
    print(f"[auto_params] markdup_memory : {params['markdup_memory']}")
    print(f"[auto_params] saved to       → {args.output}")


if __name__ == "__main__":
    main()
