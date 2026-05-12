#!/usr/bin/env python3
"""Analyze Mutect2 filtered VCF: filter breakdown and VAF distribution.

Usage:
    python scripts/analyze_vcf_filters.py \
        --vcf results/calling/WES/WES.filtered.vcf.gz \
        --bcftools-sif containers/bcftools_1.17.sif \
        [--out-prefix results/calling/WES/vcf_analysis]
"""

import argparse
import subprocess
from collections import Counter, defaultdict


def run_bcftools(args_list: list[str], sif: str) -> str:
    cmd = ["apptainer", "exec", sif, "bcftools"] + args_list
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return str(result.stdout)


def parse_vaf(format_str: str, sample_str: str) -> float:
    """Extract VAF from Mutect2 FORMAT/AF field."""
    keys = format_str.split(":")
    vals = sample_str.split(":")
    fmt = dict(zip(keys, vals, strict=False))
    if "AF" in fmt:
        try:
            return float(fmt["AF"].split(",")[0])
        except ValueError:
            pass
    if "AD" in fmt:
        try:
            ad = [int(x) for x in fmt["AD"].split(",")]
            total = sum(ad)
            return ad[1] / total if total > 0 else 0.0
        except (ValueError, IndexError):
            pass
    return -1.0


def analyze(vcf: str, sif: str) -> tuple[Counter[str], "defaultdict[str, list[float]]"]:
    # Get all variants (including filtered) with FILTER and FORMAT
    cmd = ["view", "-H", vcf]
    output = run_bcftools(cmd, sif)

    filter_counts: Counter[str] = Counter()
    filter_vafs: defaultdict[str, list[float]] = defaultdict(list)

    for line in output.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 10:
            continue
        _chrom, _pos, _id_, ref, alt, _qual, filt, _info, fmt = fields[:9]
        sample = fields[9]

        # Skip non-SNV
        if len(ref) != 1 or len(alt) != 1 or alt == ".":
            continue

        filters = filt.split(";")
        vaf = parse_vaf(fmt, sample)

        for f in filters:
            filter_counts[f] += 1
            if vaf >= 0:
                filter_vafs[f].append(vaf)

    return filter_counts, filter_vafs


def vaf_histogram(vafs: list[float], bins: int = 20) -> str:
    if not vafs:
        return "  (no data)"
    min_v, max_v = 0.0, 1.0
    step = (max_v - min_v) / bins
    counts = [0] * bins
    for v in vafs:
        idx = min(int(v / step), bins - 1)
        counts[idx] += 1
    max_count = max(counts) if counts else 1
    lines = []
    for i, c in enumerate(counts):
        lo = i * step
        hi = lo + step
        bar = "#" * int(c / max_count * 40)
        lines.append(f"  {lo:.2f}-{hi:.2f} | {bar:<40} {c}")
    return "\n".join(lines)


def percentiles(vafs: list[float]) -> str:
    if not vafs:
        return "N/A"
    s = sorted(vafs)
    n = len(s)

    def p(q: float) -> float:
        return s[int(n * q)]

    return (
        f"min={s[0]:.4f} p10={p(0.1):.4f} p25={p(0.25):.4f} "
        f"median={p(0.5):.4f} p75={p(0.75):.4f} p90={p(0.9):.4f} max={s[-1]:.4f}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--bcftools-sif", required=True, dest="bcftools_sif")
    parser.add_argument("--out-prefix", default=None, dest="out_prefix")
    args = parser.parse_args()

    print("=" * 60)
    print("VCF Filter Analysis")
    print(f"VCF: {args.vcf}")
    print("=" * 60)

    filter_counts, filter_vafs = analyze(args.vcf, args.bcftools_sif)

    total_snv = sum(filter_counts.values())
    pass_count = filter_counts.get("PASS", 0)

    print(f"\n[Filter breakdown] total SNVs = {total_snv}")
    print(f"{'FILTER':<45} {'COUNT':>8}  {'%':>6}")
    print("-" * 62)
    for f, c in sorted(filter_counts.items(), key=lambda x: -x[1]):
        pct = c / total_snv * 100 if total_snv else 0
        print(f"  {f:<43} {c:>8}  {pct:>5.1f}%")

    print("\n[VAF distribution by filter status]")
    for filt in ["PASS"] + [f for f in sorted(filter_counts) if f != "PASS"]:
        vafs = filter_vafs.get(filt, [])
        if not vafs:
            continue
        print(f"\n--- {filt} (n={len(vafs)}) ---")
        print(f"  {percentiles(vafs)}")
        print(vaf_histogram(vafs))

    # Save TSV for downstream use
    if args.out_prefix:
        tsv_path = args.out_prefix + ".filter_stats.tsv"
        with open(tsv_path, "w") as fh:
            fh.write("filter\tcount\tpct\tvaf_median\tvaf_mean\n")
            for filt, c in sorted(filter_counts.items(), key=lambda x: -x[1]):
                vafs = filter_vafs.get(filt, [])
                pct = c / total_snv * 100 if total_snv else 0
                median = sorted(vafs)[len(vafs) // 2] if vafs else -1.0
                mean = sum(vafs) / len(vafs) if vafs else -1.0
                fh.write(f"{filt}\t{c}\t{pct:.2f}\t{median:.4f}\t{mean:.4f}\n")
        print(f"\n[Saved] {tsv_path}")

    print("\n" + "=" * 60)
    print(f"PASS SNVs: {pass_count} / {total_snv} ({pass_count / total_snv * 100:.2f}%)")
    print("=" * 60)


if __name__ == "__main__":
    main()
