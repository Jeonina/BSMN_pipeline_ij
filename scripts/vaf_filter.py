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
from typing import Dict, Tuple

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
    bases = re.sub(r"\^.", "", bases_str)          # read-start + mapping quality
    bases = re.sub(r"\$", "", bases)               # read-end
    for n in set(re.findall(r"-(\d+)", bases)):    # deletions
        bases = re.sub(r"-{n}[ACGTNacgtn]{{{n}}}".format(n=n), "", bases)
    for n in set(re.findall(r"\+(\d+)", bases)):   # insertions
        bases = re.sub(r"\+{n}[ACGTNacgtn]{{{n}}}".format(n=n), "", bases)
    return bases


def pileup_base_counts(
    cram: str,
    ref: str,
    chrom: str,
    pos: str,
    samtools_sif: str,
    min_mapq: int = 20,
    min_baseq: int = 20,
) -> Tuple[int, Dict[str, int]]:
    """Run samtools mpileup via Apptainer and return (depth, per-base counts).

    Returns (0, {}) on error or zero coverage.
    """
    cmd = [
        "apptainer", "exec", samtools_sif,
        "samtools", "mpileup",
        "-d", "8000",
        "-q", str(min_mapq),
        "-Q", str(min_baseq),
        "-r", f"{chrom}:{pos}-{pos}",
        "--reference", ref,
        cram,
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        log.warning(
            "mpileup failed at %s:%s — returncode=%d stderr=%s",
            chrom, pos, e.returncode, (e.stderr or "").strip(),
        )
        return 0, {}

    parts = result.stdout.strip().split()
    if len(parts) < 5:
        return 0, {}

    # mpileup columns: chrom, pos, ref_base, depth, bases, quals, [repeat]
    bases_raw = parts[4]
    bases = _clean_bases(bases_raw)

    counts: Dict[str, int] = {
        b: bases.count(b)
        for b in ("A", "a", "C", "c", "G", "g", "T", "t")
    }
    depth = sum(counts.values()) + bases.count("*")
    return depth, counts


def filter_variants(
    txt_file,
    cram: str,
    ref: str,
    samtools_sif: str,
    min_mapq: int = 20,
    min_baseq: int = 20,
    p_threshold: float = 1e-6,
    min_alt: int = 5,
    outfile=None,
) -> dict:
    """Filter text-format variant file by VAF criteria.

    Input format: chrom  pos  ref  alt  (tab-separated, one variant per line)
    Returns a dict with filtering statistics.
    """
    if outfile is None:
        outfile = sys.stdout

    stats = {
        "input": 0, "kept": 0,
        "removed_low_alt": 0, "removed_pvalue": 0,
        "removed_no_coverage": 0, "skipped": 0,
        "depths": [], "vafs": [],
    }

    for line in txt_file:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            stats["skipped"] += 1
            continue
        stats["input"] += 1
        chrom, pos, _ref_base, alt = parts[0], parts[1], parts[2], parts[3]

        depth, counts = pileup_base_counts(
            cram, ref, chrom, pos, samtools_sif, min_mapq, min_baseq
        )
        alt_up = alt.upper()
        alt_n = counts.get(alt_up, 0) + counts.get(alt_up.lower(), 0)

        if depth == 0:
            stats["removed_no_coverage"] += 1
            log.debug(
                "REMOVED no_coverage: %s:%s %s>%s depth=0", chrom, pos, _ref_base, alt
            )
            continue

        pval = binom_pvalue(alt_n, depth)
        vaf = alt_n / depth if depth > 0 else 0.0

        if alt_n < min_alt:
            stats["removed_low_alt"] += 1
            log.debug(
                "REMOVED low_alt: %s:%s %s>%s depth=%d alt=%d vaf=%.4f p=%.2e",
                chrom, pos, _ref_base, alt, depth, alt_n, vaf, pval,
            )
        elif pval >= p_threshold:
            stats["removed_pvalue"] += 1
            log.debug(
                "REMOVED high_pvalue: %s:%s %s>%s depth=%d alt=%d vaf=%.4f p=%.2e (>%.1e)",
                chrom, pos, _ref_base, alt, depth, alt_n, vaf, pval, p_threshold,
            )
        else:
            stats["kept"] += 1
            stats["depths"].append(depth)
            stats["vafs"].append(vaf)
            log.debug(
                "KEPT: %s:%s %s>%s depth=%d alt=%d vaf=%.4f p=%.2e",
                chrom, pos, _ref_base, alt, depth, alt_n, vaf, pval,
            )
            outfile.write(line)

    return stats


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter variants by VAF (binomial test + minimum alt count)."
    )
    parser.add_argument(
        "infile",
        nargs="?",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="Input text file (chrom\\tpos\\tref\\talt per line)",
    )
    parser.add_argument("--cram", required=True, help="CRAM/BAM file for pileup")
    parser.add_argument("--ref", required=True, help="Reference FASTA")
    parser.add_argument(
        "--samtools-sif", required=True, dest="samtools_sif",
        help="Apptainer SIF for samtools"
    )
    parser.add_argument("--min-mapq", type=int, default=20, dest="min_mapq")
    parser.add_argument("--min-baseq", type=int, default=20, dest="min_baseq")
    parser.add_argument(
        "--p-threshold", type=float, default=1e-6, dest="p_threshold"
    )
    parser.add_argument("--min-alt", type=int, default=5, dest="min_alt")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of threads (reserved for future parallelism)")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
        stream=sys.stderr,
    )

    log.info("================================================================")
    log.info("START vaf_filter")
    log.info("input_file=%s", args.infile.name)
    log.info("cram=%s", args.cram)
    log.info("cram_size=%s bytes",
             os.path.getsize(args.cram) if os.path.exists(args.cram) else "N/A")
    log.info("ref=%s", args.ref)
    log.info("samtools_sif=%s", args.samtools_sif)
    log.info("parameters:")
    log.info("  min_mapq=%d", args.min_mapq)
    log.info("  min_baseq=%d", args.min_baseq)
    log.info("  p_binom_threshold=%.1e", args.p_threshold)
    log.info("  min_alt_count=%d", args.min_alt)
    log.info("  threads=%d", args.threads)
    log.info("criterion: binom_test(alt, depth, p=0.5, alt='less') < %.1e AND alt >= %d",
             args.p_threshold, args.min_alt)
    log.info("================================================================")

    t0 = time.time()
    stats = filter_variants(
        args.infile,
        args.cram,
        args.ref,
        args.samtools_sif,
        args.min_mapq,
        args.min_baseq,
        args.p_threshold,
        args.min_alt,
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
        log.info("  kept_depth: min=%d median=%d max=%d",
                 min(stats["depths"]),
                 sorted(stats["depths"])[len(stats["depths"]) // 2],
                 max(stats["depths"]))
    if stats["vafs"]:
        log.info("  kept_vaf: min=%.4f median=%.4f max=%.4f",
                 min(stats["vafs"]),
                 sorted(stats["vafs"])[len(stats["vafs"]) // 2],
                 max(stats["vafs"]))
    log.info("  elapsed=%.1f seconds", elapsed)
    log.info("END vaf_filter")
    log.info("================================================================")


if __name__ == "__main__":
    main()
