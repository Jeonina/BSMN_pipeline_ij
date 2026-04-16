#!/usr/bin/env python3
"""Accessibility filter using 1KG strict mask FASTA.

Extracts PASS SNVs from a Mutect2-filtered VCF and keeps only positions
where the 1KG strict mask base is 'P' (accessible).

Uses samtools faidx via Apptainer for mask queries, parallelized with
ThreadPoolExecutor.

Usage:
    python scripts/accessibility_filter.py \
        --vcf sample.filtered.vcf.gz \
        --mask 1KG.20160622.strict_mask.hg38_GRCh38.fa.gz \
        --bcftools-sif containers/bcftools.sif \
        --samtools-sif containers/samtools.sif \
        --threads 20 \
        > output.txt
"""

import argparse
import logging
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from typing import List, Tuple

log = logging.getLogger("accessibility_filter")


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


def query_mask(args_tuple: Tuple) -> Tuple[Tuple, str]:
    """Query mask FASTA at chrom:pos. Returns ((chrom,pos,ref,alt), mask_base)."""
    chrom, pos, ref, alt, mask, samtools_sif = args_tuple
    region = f"{chrom}:{pos}-{pos}"
    cmd = [
        "apptainer", "exec", samtools_sif,
        "samtools", "faidx", mask, region,
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        lines = result.stdout.strip().splitlines()
        base = lines[-1].strip() if len(lines) >= 2 else "?"
    except subprocess.CalledProcessError as e:
        log.warning(
            "faidx failed at %s:%s — returncode=%d stderr=%s",
            chrom, pos, e.returncode, (e.stderr or "").strip(),
        )
        base = "?"
    return (chrom, pos, ref, alt), base


def filter_variants(
    vcf: str,
    mask: str,
    bcftools_sif: str,
    samtools_sif: str,
    outfile=None,
    n_threads: int = 1,
) -> dict:
    if outfile is None:
        outfile = sys.stdout

    stats = {"input": 0, "kept": 0, "removed": 0, "mask_errors": 0}

    log.info("Extracting PASS SNVs from VCF...")
    variants = extract_pass_snvs(vcf, bcftools_sif)
    stats["input"] = len(variants)
    log.info("  total_PASS_SNVs=%d", stats["input"])

    # Build args tuples for executor.map
    work_items = [(chrom, pos, ref, alt, mask, samtools_sif) for chrom, pos, ref, alt in variants]

    log.info("Querying mask with %d threads...", n_threads)
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        results = list(executor.map(query_mask, work_items))

    for (chrom, pos, ref, alt), mask_base in results:
        if mask_base == "?":
            stats["mask_errors"] += 1
        if mask_base == "P":
            stats["kept"] += 1
            log.debug("KEPT: %s:%s %s>%s mask=%s", chrom, pos, ref, alt, mask_base)
            outfile.write(f"{chrom}\t{pos}\t{ref}\t{alt}\n")
        else:
            stats["removed"] += 1
            log.debug("REMOVED: %s:%s %s>%s mask=%s (not 'P')", chrom, pos, ref, alt, mask_base)

    return stats


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Accessibility filter: keep PASS SNVs where 1KG mask base == 'P'."
    )
    parser.add_argument("--vcf", required=True, help="Mutect2-filtered VCF (gzipped + tabix)")
    parser.add_argument("--mask", required=True,
                        help="1KG strict mask FASTA (.fa.gz, faidx indexed)")
    parser.add_argument("--bcftools-sif", required=True, dest="bcftools_sif",
                        help="Apptainer SIF for bcftools")
    parser.add_argument("--samtools-sif", required=True, dest="samtools_sif",
                        help="Apptainer SIF for samtools")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of parallel threads for faidx queries")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
        stream=sys.stderr,
    )

    log.info("================================================================")
    log.info("START accessibility_filter")
    log.info("vcf=%s", args.vcf)
    log.info("mask=%s", args.mask)
    log.info("bcftools_sif=%s", args.bcftools_sif)
    log.info("samtools_sif=%s", args.samtools_sif)
    log.info("threads=%d", args.threads)
    log.info("criterion: 1KG mask base == 'P' (accessible)")
    log.info("================================================================")

    t0 = time.time()
    stats = filter_variants(
        args.vcf, args.mask, args.bcftools_sif, args.samtools_sif,
        n_threads=args.threads,
    )
    sys.stdout.flush()
    elapsed = time.time() - t0

    log.info("================================================================")
    log.info("RESULTS:")
    log.info("  variants_input=%d", stats["input"])
    log.info("  variants_kept=%d", stats["kept"])
    log.info("  variants_removed=%d", stats["removed"])
    log.info("  mask_errors=%d", stats["mask_errors"])
    if stats["input"] > 0:
        log.info("  pass_rate=%.1f%%", stats["kept"] / stats["input"] * 100)
    log.info("  elapsed=%.1f seconds", elapsed)
    log.info("END accessibility_filter")
    log.info("================================================================")


if __name__ == "__main__":
    main()
