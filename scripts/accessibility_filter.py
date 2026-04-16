#!/usr/bin/env python3
"""Accessibility filter using 1KG strict mask BED.

Extracts PASS SNVs from a Mutect2-filtered VCF and keeps only positions
that fall within accessible regions defined by the 1KG strict mask BED file.

BED intervals are 0-based half-open [start, end). VCF positions are 1-based.
A variant at VCF position P is accessible if any BED interval satisfies:
  start < P <= end  (i.e., start+1 <= P <= end in 1-based coordinates)

Uses bcftools via Apptainer for VCF parsing. No samtools dependency.

Usage:
    python scripts/accessibility_filter.py \
        --vcf sample.filtered.vcf.gz \
        --bed 1KG.20160622.strict_mask.hg38_GRCh38.bed \
        --bcftools-sif containers/bcftools.sif \
        > output.txt
"""

import argparse
import bisect
import logging
import subprocess
import sys
import time
from collections import defaultdict
from typing import Dict, List, Tuple

log = logging.getLogger("accessibility_filter")


def load_bed(bed_path: str) -> Dict[str, Tuple[List[int], List[int]]]:
    """Load BED file into memory as sorted interval arrays per chromosome.

    Returns dict mapping chrom -> (starts, ends) where both lists are
    sorted by start position for binary search.
    """
    intervals: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            intervals[chrom].append((start, end))

    # Sort by start and split into parallel arrays for bisect
    result = {}
    for chrom, ivs in intervals.items():
        ivs.sort()
        starts = [s for s, _ in ivs]
        ends = [e for _, e in ivs]
        result[chrom] = (starts, ends)

    return result


def is_accessible(chrom: str, pos_1based: int,
                  bed: Dict[str, Tuple[List[int], List[int]]]) -> bool:
    """Return True if pos_1based falls within any BED interval on chrom.

    BED is 0-based half-open [start, end). VCF pos is 1-based.
    Accessible when: start < pos_1based <= end
    """
    if chrom not in bed:
        return False
    starts, ends = bed[chrom]
    pos_0based = pos_1based - 1  # convert to 0-based

    # Find the rightmost interval whose start <= pos_0based
    idx = bisect.bisect_right(starts, pos_0based) - 1
    if idx < 0:
        return False
    return pos_0based < ends[idx]


def extract_pass_snvs(vcf: str, bcftools_sif: str) -> List[Tuple[str, str, str, str]]:
    """Extract PASS SNVs from VCF using bcftools. Returns list of (chrom, pos, ref, alt)."""
    cmd = [
        "apptainer", "exec", bcftools_sif,
        "bcftools", "view", "-H", "-f", "PASS", "-v", "snps", vcf,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    variants = []
    for line in result.stdout.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
        if "," in alt:
            continue
        variants.append((chrom, pos, ref, alt))
    return variants


def filter_variants(
    vcf: str,
    bed_path: str,
    bcftools_sif: str,
    outfile=None,
) -> dict:
    if outfile is None:
        outfile = sys.stdout

    stats = {"input": 0, "kept": 0, "removed": 0, "no_chrom": 0}

    log.info("Loading BED mask into memory...")
    bed = load_bed(bed_path)
    log.info("  chromosomes_in_bed=%d", len(bed))

    log.info("Extracting PASS SNVs from VCF...")
    variants = extract_pass_snvs(vcf, bcftools_sif)
    stats["input"] = len(variants)
    log.info("  total_PASS_SNVs=%d", stats["input"])

    log.info("Filtering variants against BED mask...")
    for chrom, pos, ref, alt in variants:
        if chrom not in bed:
            stats["no_chrom"] += 1
            stats["removed"] += 1
            log.debug("REMOVED: %s:%s %s>%s — chrom not in BED", chrom, pos, ref, alt)
            continue
        if is_accessible(chrom, int(pos), bed):
            stats["kept"] += 1
            log.debug("KEPT: %s:%s %s>%s", chrom, pos, ref, alt)
            outfile.write(f"{chrom}\t{pos}\t{ref}\t{alt}\n")
        else:
            stats["removed"] += 1
            log.debug("REMOVED: %s:%s %s>%s — not in accessible region", chrom, pos, ref, alt)

    return stats


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Accessibility filter: keep PASS SNVs within 1KG strict mask BED regions."
    )
    parser.add_argument("--vcf", required=True, help="Mutect2-filtered VCF (gzipped + tabix)")
    parser.add_argument("--bed", required=True,
                        help="1KG strict mask BED file (0-based half-open intervals)")
    parser.add_argument("--bcftools-sif", required=True, dest="bcftools_sif",
                        help="Apptainer SIF for bcftools")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
        stream=sys.stderr,
    )

    log.info("================================================================")
    log.info("START accessibility_filter")
    log.info("vcf=%s", args.vcf)
    log.info("bed=%s", args.bed)
    log.info("bcftools_sif=%s", args.bcftools_sif)
    log.info("criterion: variant position within 1KG strict mask BED interval")
    log.info("================================================================")

    t0 = time.time()
    stats = filter_variants(args.vcf, args.bed, args.bcftools_sif)
    sys.stdout.flush()
    elapsed = time.time() - t0

    log.info("================================================================")
    log.info("RESULTS:")
    log.info("  variants_input=%d", stats["input"])
    log.info("  variants_kept=%d", stats["kept"])
    log.info("  variants_removed=%d", stats["removed"])
    log.info("  chrom_not_in_bed=%d", stats["no_chrom"])
    if stats["input"] > 0:
        log.info("  pass_rate=%.1f%%", stats["kept"] / stats["input"] * 100)
    log.info("  elapsed=%.1f seconds", elapsed)
    log.info("END accessibility_filter")
    log.info("================================================================")


if __name__ == "__main__":
    main()
