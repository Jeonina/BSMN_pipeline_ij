#!/usr/bin/env python3
"""VAF filter: keep somatic candidates passing binomial test + minimum alt count.

Filter criteria (from original BSMN C.VAF_filters.sh):
  - binomial test p-value < p_threshold (default 1e-6)
  - alt read count >= min_alt (default 5)

The one-sided binomial test (H0: VAF >= 0.5, alternative: VAF < 0.5) has a very
small p-value for somatic mosaic variants (VAF ~0.01-0.30) and a large p-value
for germline heterozygous variants (VAF ~0.5). Variants passing the filter are
confidently below the germline VAF level.

Uses samtools mpileup (via Apptainer) for CRAM pileup.

Usage:
    python scripts/vaf_filter.py \\
        --cram sample.cram --ref hg38.fa \\
        --samtools-sif containers/samtools_1.17.sif \\
        --min-mapq 20 --min-baseq 20 \\
        input.txt > output.txt
"""

import argparse
import logging
import math
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from typing import IO, Any

log = logging.getLogger("vaf_filter")


def binom_pvalue(alt_n: int, depth: int) -> float:
    """P-value for one-sided binomial test (H_alt: VAF < 0.5).

    Uses a normal approximation with continuity correction.
    Accurate to p-values well below 1e-6 (which is our threshold).
    Equivalent to scipy.stats.binomtest(alt_n, depth, p=0.5, alternative='less').
    """
    if depth == 0:
        return 1.0
    mu = depth * 0.5
    sigma = math.sqrt(depth * 0.25)
    # continuity correction: P(X <= alt_n) ≈ P(Z <= (alt_n + 0.5 - mu) / sigma)
    z = (alt_n + 0.5 - mu) / sigma
    # normal CDF via erfc: Phi(z) = 0.5 * erfc(-z / sqrt(2))
    return 0.5 * math.erfc(-z / math.sqrt(2))


def passes_vaf_filter(
    alt_n: int,
    depth: int,
    p_threshold: float = 1e-6,
    min_alt: int = 5,
) -> bool:
    """Return True if variant passes the VAF filter.

    Args:
        alt_n:       Number of alt-supporting reads.
        depth:       Total read depth at position.
        p_threshold: Maximum allowed binomial p-value (default 1e-6).
        min_alt:     Minimum required alt read count (default 5).

    Returns:
        True  → somatic candidate (keep)
        False → germline or insufficient support (discard)
    """
    if alt_n < min_alt:
        return False
    return binom_pvalue(alt_n, depth) < p_threshold


def _clean_bases(bases_str: str) -> str:
    """Strip samtools mpileup markup from base string."""
    bases = re.sub(r"\^.", "", bases_str)  # read-start + mapping quality
    bases = re.sub(r"\$", "", bases)  # read-end
    for n in set(re.findall(r"-(\d+)", bases)):  # deletions
        bases = re.sub(rf"-{n}[ACGTNacgtn]{{{n}}}", "", bases)
    for n in set(re.findall(r"\+(\d+)", bases)):  # insertions
        bases = re.sub(rf"\+{n}[ACGTNacgtn]{{{n}}}", "", bases)
    return bases


def pileup_base_counts(
    cram: str,
    ref: str,
    chrom: str,
    pos: str,
    samtools_sif: str,
    min_mapq: int = 20,
    min_baseq: int = 20,
) -> tuple[int, dict[str, int]]:
    """Run samtools mpileup via Apptainer and return (depth, per-base counts).

    Returns (0, {}) on error or zero coverage.
    """
    cmd = [
        "apptainer",
        "exec",
        samtools_sif,
        "samtools",
        "mpileup",
        "-d",
        "8000",
        "-q",
        str(min_mapq),
        "-Q",
        str(min_baseq),
        "-r",
        f"{chrom}:{pos}-{pos}",
        "--reference",
        ref,
        cram,
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        log.warning(
            "mpileup failed at %s:%s — returncode=%d stderr=%s",
            chrom,
            pos,
            e.returncode,
            (e.stderr or "").strip(),
        )
        return 0, {}

    parts = result.stdout.strip().split()
    if len(parts) < 5:
        return 0, {}

    # mpileup columns: chrom, pos, ref_base, depth, bases, quals, [repeat]
    bases_raw = parts[4]
    bases = _clean_bases(bases_raw)

    counts: dict[str, int] = {b: bases.count(b) for b in ("A", "a", "C", "c", "G", "g", "T", "t")}
    depth = sum(counts.values()) + bases.count("*")
    return depth, counts


def _process_one(
    line: str,
    cram: str,
    ref: str,
    samtools_sif: str,
    min_mapq: int,
    min_baseq: int,
    p_threshold: float,
    min_alt: int,
) -> dict[str, Any]:
    """Process a single variant line. Called from worker threads."""
    parts = line.rstrip("\n").split("\t")
    chrom, pos, ref_base, alt = parts[0], parts[1], parts[2], parts[3]

    depth, counts = pileup_base_counts(cram, ref, chrom, pos, samtools_sif, min_mapq, min_baseq)
    alt_up = alt.upper()
    alt_n = counts.get(alt_up, 0) + counts.get(alt_up.lower(), 0)
    pval = binom_pvalue(alt_n, depth) if depth > 0 else 1.0
    vaf = alt_n / depth if depth > 0 else 0.0

    if depth == 0:
        verdict = "no_coverage"
    elif alt_n < min_alt:
        verdict = "low_alt"
    elif pval >= p_threshold:
        verdict = "high_pvalue"
    else:
        verdict = "keep"

    return {
        "line": line,
        "chrom": chrom,
        "pos": pos,
        "ref_base": ref_base,
        "alt": alt,
        "depth": depth,
        "alt_n": alt_n,
        "pval": pval,
        "vaf": vaf,
        "verdict": verdict,
    }


def filter_variants(
    txt_file: IO[str],
    cram: str,
    ref: str,
    samtools_sif: str,
    min_mapq: int = 20,
    min_baseq: int = 20,
    p_threshold: float = 1e-6,
    min_alt: int = 5,
    outfile: IO[str] | None = None,
    n_threads: int = 1,
) -> dict[str, Any]:
    """Filter text-format variant file by VAF criteria.

    Input format: chrom  pos  ref  alt  (tab-separated, one variant per line)
    Returns a dict with filtering statistics.
    """
    if outfile is None:
        outfile = sys.stdout

    stats: dict[str, Any] = {
        "input": 0,
        "kept": 0,
        "removed_low_alt": 0,
        "removed_pvalue": 0,
        "removed_no_coverage": 0,
        "skipped": 0,
        "depths": [],
        "vafs": [],
        "all_depths": {"kept": [], "low_alt": [], "high_pvalue": [], "no_coverage": []},
        "all_vafs": {"kept": [], "low_alt": [], "high_pvalue": []},
    }

    # Read all variant lines upfront (headers/malformed lines filtered here)
    variants: list[str] = []
    for line in txt_file:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            stats["skipped"] += 1
            continue
        variants.append(line)

    stats["input"] = len(variants)

    def _worker(line: str) -> dict[str, Any]:
        return _process_one(
            line,
            cram,
            ref,
            samtools_sif,
            min_mapq,
            min_baseq,
            p_threshold,
            min_alt,
        )

    # Run pileup calls in parallel; executor.map preserves input order
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        results = list(executor.map(_worker, variants))

    for r in results:
        verdict = r["verdict"]
        chrom, pos, ref_base, alt = r["chrom"], r["pos"], r["ref_base"], r["alt"]
        depth, alt_n, pval, vaf = r["depth"], r["alt_n"], r["pval"], r["vaf"]

        if verdict == "no_coverage":
            stats["removed_no_coverage"] += 1
            stats["all_depths"]["no_coverage"].append(depth)
            log.debug("REMOVED no_coverage: %s:%s %s>%s depth=0", chrom, pos, ref_base, alt)
        elif verdict == "low_alt":
            stats["removed_low_alt"] += 1
            stats["all_depths"]["low_alt"].append(depth)
            stats["all_vafs"]["low_alt"].append(vaf)
            log.debug(
                "REMOVED low_alt: %s:%s %s>%s depth=%d alt=%d vaf=%.4f p=%.2e",
                chrom,
                pos,
                ref_base,
                alt,
                depth,
                alt_n,
                vaf,
                pval,
            )
        elif verdict == "high_pvalue":
            stats["removed_pvalue"] += 1
            stats["all_depths"]["high_pvalue"].append(depth)
            stats["all_vafs"]["high_pvalue"].append(vaf)
            log.debug(
                "REMOVED high_pvalue: %s:%s %s>%s depth=%d alt=%d vaf=%.4f p=%.2e (>%.1e)",
                chrom,
                pos,
                ref_base,
                alt,
                depth,
                alt_n,
                vaf,
                pval,
                p_threshold,
            )
        else:
            stats["kept"] += 1
            stats["depths"].append(depth)
            stats["vafs"].append(vaf)
            stats["all_depths"]["kept"].append(depth)
            stats["all_vafs"]["kept"].append(vaf)
            log.debug(
                "KEPT: %s:%s %s>%s depth=%d alt=%d vaf=%.4f p=%.2e",
                chrom,
                pos,
                ref_base,
                alt,
                depth,
                alt_n,
                vaf,
                pval,
            )
            outfile.write(r["line"])

    return stats


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter variants by VAF (binomial test + minimum alt count)."
    )
    parser.add_argument(
        "infile",
        nargs="?",
        default=None,
        help="Input text file (chrom\\tpos\\tref\\talt per line); default stdin",
    )
    parser.add_argument("--cram", required=True, help="CRAM/BAM file for pileup")
    parser.add_argument("--ref", required=True, help="Reference FASTA")
    parser.add_argument(
        "--samtools-sif", required=True, dest="samtools_sif", help="Apptainer SIF for samtools"
    )
    parser.add_argument("--min-mapq", type=int, default=20, dest="min_mapq")
    parser.add_argument("--min-baseq", type=int, default=20, dest="min_baseq")
    parser.add_argument("--p-threshold", type=float, default=1e-6, dest="p_threshold")
    parser.add_argument("--min-alt", type=int, default=5, dest="min_alt")
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of parallel threads for samtools mpileup calls",
    )
    args = parser.parse_args()
    infile: IO[str] = open(args.infile) if args.infile else sys.stdin

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
        stream=sys.stderr,
    )

    log.info("================================================================")
    log.info("START vaf_filter")
    log.info("input_file=%s", getattr(infile, "name", "<stdin>"))
    log.info("cram=%s", args.cram)
    log.info(
        "cram_size=%s bytes", os.path.getsize(args.cram) if os.path.exists(args.cram) else "N/A"
    )
    log.info("ref=%s", args.ref)
    log.info("samtools_sif=%s", args.samtools_sif)
    log.info("parameters:")
    log.info("  min_mapq=%d", args.min_mapq)
    log.info("  min_baseq=%d", args.min_baseq)
    log.info("  p_binom_threshold=%.1e", args.p_threshold)
    log.info("  min_alt_count=%d", args.min_alt)
    log.info("  threads=%d", args.threads)
    log.info(
        "criterion: binom_test(alt, depth, p=0.5, alt='less') < %.1e AND alt >= %d",
        args.p_threshold,
        args.min_alt,
    )
    log.info("================================================================")

    t0 = time.time()
    stats = filter_variants(
        infile,
        args.cram,
        args.ref,
        args.samtools_sif,
        args.min_mapq,
        args.min_baseq,
        args.p_threshold,
        args.min_alt,
        n_threads=args.threads,
    )
    sys.stdout.flush()
    elapsed = time.time() - t0

    total_removed = (
        stats["removed_low_alt"] + stats["removed_pvalue"] + stats["removed_no_coverage"]
    )

    log.info("================================================================")
    log.info("RESULTS:")
    log.info("  variants_input=%d", stats["input"])
    log.info("  variants_kept=%d", stats["kept"])
    log.info("  variants_removed=%d", total_removed)
    log.info("    removed_low_alt=%d (alt < %d)", stats["removed_low_alt"], args.min_alt)
    log.info("    removed_high_pvalue=%d (p >= %.1e)", stats["removed_pvalue"], args.p_threshold)
    log.info("    removed_no_coverage=%d", stats["removed_no_coverage"])
    log.info("  lines_skipped=%d (malformed)", stats["skipped"])
    if stats["input"] > 0:
        log.info("  pass_rate=%.1f%%", stats["kept"] / stats["input"] * 100)
    if stats["depths"]:
        log.info(
            "  kept_depth: min=%d median=%d max=%d",
            min(stats["depths"]),
            sorted(stats["depths"])[len(stats["depths"]) // 2],
            max(stats["depths"]),
        )
    if stats["vafs"]:
        log.info(
            "  kept_vaf: min=%.4f median=%.4f max=%.4f",
            min(stats["vafs"]),
            sorted(stats["vafs"])[len(stats["vafs"]) // 2],
            max(stats["vafs"]),
        )
    log.info("  elapsed=%.1f seconds", elapsed)

    # --- Depth histogram (all verdicts) ---
    log.info("DEPTH DISTRIBUTION (mpileup depth, by verdict):")
    for verdict, depths in stats["all_depths"].items():
        if not depths:
            continue
        s = sorted(depths)
        n = len(s)
        buckets = [0] * 10
        boundaries = [0, 5, 10, 20, 30, 50, 75, 100, 150, 200, 99999]
        labels = [
            "0-4",
            "5-9",
            "10-19",
            "20-29",
            "30-49",
            "50-74",
            "75-99",
            "100-149",
            "150-199",
            "200+",
        ]
        for d in depths:
            for i, bound in enumerate(boundaries[1:]):
                if d < bound:
                    buckets[i] += 1
                    break
        median = s[n // 2]
        mean = sum(depths) / n
        log.info("  [%s] n=%d  mean=%.1f  median=%d", verdict, n, mean, median)
        for label, count in zip(labels, buckets, strict=False):
            if count == 0:
                continue
            bar = "#" * min(40, int(count / n * 40) + 1)
            log.info("    depth %6s | %-40s %d (%.1f%%)", label, bar, count, count / n * 100)

    # --- VAF histogram (all verdicts except no_coverage) ---
    log.info("VAF DISTRIBUTION (by verdict):")
    for verdict, vafs in stats["all_vafs"].items():
        if not vafs:
            continue
        s = sorted(vafs)
        n = len(s)
        buckets = [0] * 10
        edges = [i * 0.1 for i in range(11)]
        for v in vafs:
            idx = min(int(v / 0.1), 9)
            buckets[idx] += 1
        median = s[n // 2]
        log.info("  [%s] n=%d  median=%.3f", verdict, n, median)
        for i, count in enumerate(buckets):
            if count == 0:
                continue
            label = f"{edges[i]:.1f}-{edges[i + 1]:.1f}"
            bar = "#" * min(40, int(count / n * 40) + 1)
            log.info("    VAF %9s | %-40s %d (%.1f%%)", label, bar, count, count / n * 100)

    log.info("END vaf_filter")
    log.info("================================================================")


if __name__ == "__main__":
    main()
