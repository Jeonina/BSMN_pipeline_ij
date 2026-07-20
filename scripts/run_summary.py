#!/usr/bin/env python3
"""Pipeline run summary: collect elapsed times and variant counts from all logs.

Usage:
    python scripts/run_summary.py --sample WES --log-dir logs
    python scripts/run_summary.py --sample WES  # uses logs/ by default
"""

import argparse
import os
import re
from typing import Any

# Log file locations and patterns to extract
STEPS: list[tuple[str, str, str, list[tuple[str, str]]]] = [
    # (step_name, log_path_template, elapsed_pattern, extra_patterns)
    ("bwa_map", "logs/mapping/{sample}/bwa_mem.log", r"elapsed=(\d+\.\d+) seconds", []),
    ("mark_dup", "logs/mapping/{sample}/markdup.log", r"elapsed=(\d+\.\d+) seconds", []),
    (
        "mutect2_scatter",
        "logs/calling/{sample}/mutect2_scatter.chr1.log",
        r"Elapsed time: ([\d.]+) minutes",
        [],
    ),
    (
        "learn_orientation",
        "logs/calling/{sample}/learn_read_orientation.log",
        r"Elapsed time: ([\d.]+) minutes",
        [],
    ),
    (
        "get_pileup",
        "logs/calling/{sample}/get_pileup_summaries.log",
        r"Elapsed time: ([\d.]+) minutes",
        [],
    ),
    (
        "calc_contamination",
        "logs/calling/{sample}/calculate_contamination.log",
        r"Elapsed time: ([\d.]+) minutes",
        [(r"contamination_value=([\d.E+-]+)", "contamination")],
    ),
    ("merge_vcfs", "logs/calling/{sample}/merge_vcfs.log", r"Elapsed time: ([\d.]+) minutes", []),
    (
        "filter_mutect",
        "logs/calling/{sample}/filter_mutect_calls.log",
        r"Elapsed time: ([\d.]+) minutes",
        [(r"Processed (\d+) total variants", "total_variants")],
    ),
    (
        "accessibility",
        "logs/filtering/{sample}/accessibility_filter.log",
        r"elapsed=([\d.]+) seconds",
        [
            (r"variants_input=(\d+)", "input"),
            (r"variants_kept=(\d+)", "kept"),
            (r"pass_rate=([\d.]+)%", "pass_rate"),
        ],
    ),
    (
        "germline",
        "logs/filtering/{sample}/germline_filter.log",
        r"elapsed=([\d.]+) seconds",
        [(r"variants_kept=(\d+)", "kept"), (r"pass_rate=([\d.]+)%", "pass_rate")],
    ),
    (
        "vaf",
        "logs/filtering/{sample}/vaf_filter.log",
        r"elapsed=([\d.]+) seconds",
        [
            (r"variants_input=(\d+)", "input"),
            (r"variants_kept=(\d+)", "kept"),
            (r"pass_rate=([\d.]+)%", "pass_rate"),
            (r"removed_low_alt=(\d+)", "low_alt"),
            (r"removed_high_pvalue=(\d+)", "high_pvalue"),
            (r"removed_no_coverage=(\d+)", "no_coverage"),
        ],
    ),
    (
        "cnvnator",
        "logs/filtering/{sample}/cnvnator_filter.log",
        r"elapsed=([\d.]+) seconds",
        [(r"cnvnator: (\d+) ->", "input"), (r"kept=(\d+)", "kept")],
    ),
    # mayo is the BSMN E-step alternative to MosaicForecast (two options, not a
    # chain). The default cascade uses MosaicForecast only, so mayo_filter is not
    # run and produces no log. Re-add this entry if you wire mayo back in.
    (
        "mosaicforecast",
        "logs/filtering/{sample}/mosaicforecast_filter.log",
        r"elapsed=([\d.]+) seconds",
        [(r"mosaicforecast: (\d+) ->", "input"), (r"variants_kept=(\d+)", "kept")],
    ),
    (
        "pon_mask",
        "logs/filtering/{sample}/pon_mask_filter.log",
        r"elapsed=([\d.]+) seconds",
        [(r"variants_kept=(\d+)", "kept")],
    ),
]


def parse_log(path: str, elapsed_pat: str, extra_pats: list[tuple[str, str]]) -> dict[str, Any]:
    result: dict[str, Any] = {"elapsed_sec": None, "found": False}
    if not os.path.exists(path):
        return result
    result["found"] = True
    try:
        text = open(path).read()
    except OSError:
        return result

    # elapsed
    m = re.search(elapsed_pat, text)
    if m:
        val = float(m.group(1))
        # convert minutes to seconds if pattern says "minutes"
        if "minutes" in elapsed_pat:
            val *= 60
        result["elapsed_sec"] = val

    # extra fields
    for pat, key in extra_pats:
        m = re.search(pat, text)
        if m:
            try:
                result[key] = float(m.group(1))
            except ValueError:
                result[key] = m.group(1)

    return result


def fmt_time(seconds: float | None) -> str:
    if seconds is None:
        return "N/A"
    if seconds < 60:
        return f"{seconds:.1f}s"
    if seconds < 3600:
        return f"{seconds / 60:.1f}m"
    return f"{seconds / 3600:.1f}h"


def main() -> None:
    parser = argparse.ArgumentParser(description="Pipeline run summary")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--log-dir", default="logs", dest="log_dir")
    args = parser.parse_args()

    print("=" * 70)
    print(f"Pipeline Run Summary  —  sample: {args.sample}")
    print("=" * 70)
    print(f"{'Step':<22} {'Elapsed':>8}  {'Details'}")
    print("-" * 70)

    total_sec = 0.0
    for step_name, log_tmpl, elapsed_pat, extra_pats in STEPS:
        rel = log_tmpl.format(sample=args.sample).removeprefix("logs/")
        log_path = os.path.join(args.log_dir, rel)
        # support both absolute and relative path
        if not os.path.exists(log_path):
            log_path = log_tmpl.format(sample=args.sample)

        data = parse_log(log_path, elapsed_pat, extra_pats)

        if not data["found"]:
            print(f"  {step_name:<20} {'—':>8}  (log not found)")
            continue

        elapsed_str = fmt_time(data.get("elapsed_sec"))
        if data.get("elapsed_sec"):
            total_sec += data["elapsed_sec"]

        # Build details string
        details_parts = []
        if "total_variants" in data:
            details_parts.append(f"total_called={int(data['total_variants']):,}")
        if "contamination" in data:
            details_parts.append(f"contamination={data['contamination']:.4f}")
        if "input" in data and "kept" in data:
            details_parts.append(
                f"{int(data['input']):,} → {int(data['kept']):,}"
                + (f" ({data['pass_rate']:.1f}%)" if "pass_rate" in data else "")
            )
        elif "kept" in data:
            details_parts.append(f"kept={int(data['kept']):,}")
        if "low_alt" in data:
            details_parts.append(
                f"low_alt={int(data['low_alt']):,} "
                f"high_p={int(data['high_pvalue']):,} "
                f"no_cov={int(data['no_coverage']):,}"
            )
        details = "  ".join(details_parts) if details_parts else ""
        print(f"  {step_name:<20} {elapsed_str:>8}  {details}")

    print("-" * 70)
    print(f"  {'TOTAL':<20} {fmt_time(total_sec):>8}")
    print("=" * 70)


if __name__ == "__main__":
    main()
