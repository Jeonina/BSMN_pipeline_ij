#!/usr/bin/env python3
"""Assign somatic SNVs to statistically phased germline haplotypes (Bae et al. 2022).

Bae 2022 phased each individual's germline SNPs with SHAPEIT4 (1000 Genomes
reference panel + genetic map) and then asked, read by read, whether a somatic
call sits on one haplotype ("supported") or on both ("unsupported").  The
SHAPEIT4 step is plain shell (see the guide); this script does the second part.

For each sample's final SNVs it takes the donor's phased het SNPs within
--window bp, reads the read pairs over the SNV from the CRAM, and gives each
pair a haplotype by the alleles it carries at those SNPs (all informative SNPs
must agree; a conflicting pair is dropped).  Per SNV it counts alt and ref
pairs on haplotype A and B and classifies:

  supported     alt pairs on one haplotype only, and that haplotype also has
                ref pairs (a subclone of one parental chromosome = mosaic-like)
  unsupported   alt pairs on both haplotypes, at least --min-pairs on each
                (cannot arise from one mutation = artifact-like)
  germline_like alt pairs on one haplotype only and no ref pair on it
  ambiguous     alt on both, but fewer than --min-pairs on one side
  unphased      no alt pair could be given a haplotype

Runs inside the mosaicforecast image (python3 + pysam).  Outputs (--out-dir):
  haplotype_calls.tsv    one row per SNV with the pair counts and the class
  haplotype_samples.tsv  per sample: SNVs, phased, supported, unsupported, ...

Usage (VM100, pipeline directory; phased BCFs from SHAPEIT4, one column per donor):
    APPTAINERENV_REF_PATH=$PWD/resources/hg38/ref_cache/%2s/%2s/%s \\
    apptainer exec -B /mnt/data containers/mosaicforecast.sif python3 scripts/haplotype_assign.py \\
        --sample-info ~/sample_info.tsv --results results \\
        --phased-dir /mnt/data/bsmn/shapeit/phased --ref resources/hg38/Homo_sapiens_assembly38.fasta \\
        --out-dir /mnt/data/bsmn/shapeit/assign --jobs 16
"""

import argparse
import csv
import os
import sys
from collections import Counter
from multiprocessing import Pool

CLASSES = ["supported", "unsupported", "germline_like", "ambiguous", "unphased"]


def read_final(path):
    """final.txt: CHROM POS REF ALT; SNVs only."""
    out = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[1].isdigit() and len(f[2]) == 1 and len(f[3]) == 1:
                out.append((f[0], int(f[1]), f[2], f[3]))
    return out


def pair_haplotype(votes):
    """Haplotype of a read pair from its votes at phased SNPs: 'A', 'B', or None (none / conflicting)."""
    s = set(votes)
    return s.pop() if len(s) == 1 else None


def pair_allele(bases, ref, alt):
    """SNV allele of a read pair from the base(s) its mates show: 'alt', 'ref', or None (other / mates disagree)."""
    calls = {"alt" if b == alt else "ref" if b == ref else "other" for b in bases}
    if len(calls) != 1:
        return None
    c = calls.pop()
    return None if c == "other" else c


def classify(n, min_pairs):
    """n: Counter of (allele, hap) pairs -> one of CLASSES."""
    alt_a, alt_b = n[("alt", "A")], n[("alt", "B")]
    if alt_a == 0 and alt_b == 0:
        return "unphased"
    if alt_a and alt_b:
        return "unsupported" if min(alt_a, alt_b) >= min_pairs else "ambiguous"
    hap = "A" if alt_a else "B"
    return "supported" if n[("ref", hap)] > 0 else "germline_like"


def phased_snps(vf, column, chrom, start, end):
    """{pos: (allele on hap A, allele on hap B)} for the donor's phased het SNPs in [start, end] (1-based)."""
    out = {}
    for rec in vf.fetch(chrom, max(0, start - 1), end):
        if len(rec.alts or ()) != 1 or len(rec.ref) != 1 or len(rec.alts[0]) != 1:
            continue
        smp = rec.samples[column]
        gt = smp.get("GT")
        if not smp.phased or gt is None or None in gt or gt[0] == gt[1]:
            continue
        alleles = (rec.ref, rec.alts[0])
        out[rec.pos] = (alleles[gt[0]], alleles[gt[1]])
    return out


def count_pairs(bam, chrom, pos, ref, alt, snps, min_mapq, min_baseq):
    """Counter of (allele, hap) over read pairs covering the SNV (1-based pos)."""
    want = set(snps) | {pos}
    lo, hi = min(want) - 1, max(want)
    snv_bases, votes = {}, {}
    for read in bam.fetch(chrom, lo, hi):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate \
                or read.is_qcfail or read.mapping_quality < min_mapq:
            continue
        seq, qual = read.query_sequence, read.query_qualities
        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            p = rpos + 1
            if p not in want or qual[qpos] < min_baseq:
                continue
            base = seq[qpos]
            if p == pos:
                snv_bases.setdefault(read.query_name, set()).add(base)
            elif base in snps[p]:
                votes.setdefault(read.query_name, []).append("A" if base == snps[p][0] else "B")
    n = Counter()
    for name, bases in snv_bases.items():
        allele = pair_allele(bases, ref, alt)
        hap = pair_haplotype(votes.get(name, []))
        if allele and hap:
            n[(allele, hap)] += 1
    return n


def donor_column(header_samples, donor):
    cols = [s for s in header_samples if s.split("_")[0] == donor]
    return cols[0] if cols else None


def run_sample(job):
    import pysam  # inside the mosaicforecast image

    s, region, args = job
    donor = region.split("_")[0]
    calls = read_final(os.path.join(args.results, "filtering", s, f"{s}.final.txt"))
    cram = os.path.join(args.results, "mapping", s, f"{s}.cram")
    bam = pysam.AlignmentFile(cram, "rc", reference_filename=args.ref)
    rows, vfs = [], {}
    for chrom, pos, ref, alt in calls:
        path = os.path.join(args.phased_dir, f"phased_{chrom}.bcf")
        n, nsnp = Counter(), 0
        if os.path.exists(path):
            if chrom not in vfs:
                vfs[chrom] = pysam.VariantFile(path)
            col = donor_column(list(vfs[chrom].header.samples), donor)
            if col:
                snps = phased_snps(vfs[chrom], col, chrom, pos - args.window, pos + args.window)
                snps.pop(pos, None)
                nsnp = len(snps)
                if snps:
                    n = count_pairs(bam, chrom, pos, ref, alt, snps, args.min_mapq, args.min_baseq)
        rows.append([s, region, chrom, pos, ref, alt, nsnp,
                     n[("alt", "A")], n[("alt", "B")], n[("ref", "A")], n[("ref", "B")],
                     classify(n, args.min_pairs)])
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-info", required=True, help="TSV: sample_id, br_region")
    ap.add_argument("--results", default="results")
    ap.add_argument("--phased-dir", required=True, help="SHAPEIT4 output: phased_{chrom}.bcf (+ .csi)")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--window", type=int, default=1000, help="phased SNPs within this distance (bp) [1000]")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-baseq", type=int, default=20)
    ap.add_argument("--min-pairs", type=int, default=2, help="alt pairs needed on EACH haplotype for 'unsupported' [2]")
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)

    with open(args.sample_info) as fh:
        samples = [(r["sample_id"], r["br_region"]) for r in csv.DictReader(fh, delimiter="\t")]
    with Pool(args.jobs) as pool:
        per_sample = pool.map(run_sample, [(s, r, args) for s, r in samples])

    header = ["sample", "br_region", "chrom", "pos", "ref", "alt", "phased_snps",
              "alt_A", "alt_B", "ref_A", "ref_B", "class"]
    with open(os.path.join(args.out_dir, "haplotype_calls.tsv"), "w") as fo:
        fo.write("\t".join(header) + "\n")
        for rows in per_sample:
            for r in rows:
                fo.write("\t".join(map(str, r)) + "\n")
    total = Counter()
    with open(os.path.join(args.out_dir, "haplotype_samples.tsv"), "w") as fo:
        fo.write("sample\tbr_region\tsnvs\tphased\t" + "\t".join(CLASSES) + "\n")
        for (s, region), rows in zip(samples, per_sample):
            c = Counter(r[-1] for r in rows)
            total.update(c)
            fo.write(f"{s}\t{region}\t{len(rows)}\t{len(rows) - c['unphased']}\t"
                     + "\t".join(str(c[k]) for k in CLASSES) + "\n")
    print("  ".join(f"{k} {total[k]}" for k in CLASSES), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
