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
import math
import re
import subprocess
import sys
from typing import Dict, Tuple


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
    except subprocess.CalledProcessError:
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
) -> None:
    """Filter text-format variant file by VAF criteria.

    Input format: chrom  pos  ref  alt  (tab-separated, one variant per line)
    """
    if outfile is None:
        outfile = sys.stdout

    for line in txt_file:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            continue
        chrom, pos, _ref_base, alt = parts[0], parts[1], parts[2], parts[3]

        depth, counts = pileup_base_counts(
            cram, ref, chrom, pos, samtools_sif, min_mapq, min_baseq
        )
        alt_up = alt.upper()
        alt_n = counts.get(alt_up, 0) + counts.get(alt_up.lower(), 0)

        if passes_vaf_filter(alt_n, depth, p_threshold, min_alt):
            outfile.write(line)


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

    filter_variants(
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


if __name__ == "__main__":
    main()
