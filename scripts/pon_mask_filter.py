#!/usr/bin/env python3
"""PON mask filter using IUPAC FASTA.

Removes variants present in a panel-of-normals FASTA where the PON base at the
variant position matches the alt allele according to the IUPAC degenerate base
code table (from original BSMN F.PON_mask.sh).

IUPAC decision table:
  *  → always Pass (position not in PON)
  N  → always Fail
  A/C/G/T → Fail if alt matches exactly
  R  → Fail if alt ∈ {A, G}
  Y  → Fail if alt ∈ {C, T}
  S  → Fail if alt ∈ {G, C}
  W  → Fail if alt ∈ {A, T}
  K  → Fail if alt ∈ {G, T}
  M  → Fail if alt ∈ {A, C}
  B  → Fail if alt ∈ {C, G, T}
  D  → Fail if alt ∈ {A, G, T}
  H  → Fail if alt ∈ {A, C, T}
  V  → Fail if alt ∈ {A, C, G}
  default → Fail (conservative)

Usage:
    python scripts/pon_mask_filter.py \\
        --pon-fasta resources/PON.q20q20.05.5.fa \\
        --samtools-sif containers/samtools_1.17.sif \\
        input.txt > output.txt
"""

import argparse
import logging
import os
import subprocess
import sys
import time
from collections import Counter
from typing import Set

log = logging.getLogger("pon_mask_filter")

# IUPAC_TABLE maps each code to the set of nucleotides whose presence means FAIL.
# '*' maps to empty set → never matches → always Pass.
IUPAC_TABLE = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},  # N always Fail
    "*": set(),                  # * always Pass
}


def iupac_match(alt: str, iupac_base: str) -> bool:
    """Return True if alt allele matches the IUPAC code (variant is in PON → Fail).

    Case-insensitive comparison.  Unknown base codes are treated as Fail
    (conservative behaviour, same as 'default' branch in original awk).
    """
    code = iupac_base.upper()
    if code not in IUPAC_TABLE:
        return True  # unknown → conservative Fail
    return alt.upper() in IUPAC_TABLE[code]


def pon_passes(alt: str, pon_base: str) -> bool:
    """Return True if the variant is NOT masked by the PON (passes the filter)."""
    return not iupac_match(alt, pon_base)


def query_fasta(fasta: str, chrom: str, pos: str, samtools_sif: str) -> str:
    """Query a FASTA file with samtools faidx (via Apptainer).

    Returns the single base at chrom:pos-pos, or '?' on failure.
    """
    region = f"{chrom}:{pos}-{pos}"
    cmd = [
        "apptainer", "exec", samtools_sif,
        "samtools", "faidx", fasta, region,
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        lines = result.stdout.strip().splitlines()
        # faidx output: >header\nSEQUENCE
        return lines[-1].strip() if len(lines) >= 2 else "?"
    except subprocess.CalledProcessError as e:
        log.warning(
            "faidx failed at %s:%s — returncode=%d stderr=%s",
            chrom, pos, e.returncode, (e.stderr or "").strip(),
        )
        return "?"


def filter_variants(
    txt_file, pon_fasta: str, samtools_sif: str, outfile=None,
) -> dict:
    """Filter text-format variant file by PON IUPAC mask.

    Input format: chrom  pos  ref  alt  (tab-separated)
    Keeps variants where pon_passes(alt, PON_base) is True.
    Returns a dict with filtering statistics.
    """
    if outfile is None:
        outfile = sys.stdout

    stats = {
        "input": 0, "kept": 0, "removed": 0,
        "skipped": 0, "faidx_errors": 0,
    }
    pon_code_counts: Counter = Counter()

    for line in txt_file:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            stats["skipped"] += 1
            continue
        stats["input"] += 1
        chrom, pos, ref, alt = parts[0], parts[1], parts[2], parts[3]
        pon_base = query_fasta(pon_fasta, chrom, pos, samtools_sif)

        if pon_base == "?":
            stats["faidx_errors"] += 1

        pon_code_counts[pon_base.upper()] += 1

        if pon_passes(alt, pon_base):
            stats["kept"] += 1
            log.debug(
                "KEPT: %s:%s %s>%s pon_base=%s", chrom, pos, ref, alt, pon_base
            )
            outfile.write(line)
        else:
            stats["removed"] += 1
            log.debug(
                "REMOVED pon_masked: %s:%s %s>%s pon_base=%s (IUPAC match)",
                chrom, pos, ref, alt, pon_base,
            )

    stats["pon_code_distribution"] = dict(pon_code_counts.most_common())
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(
        description="PON mask filter: remove variants present in panel-of-normals FASTA."
    )
    parser.add_argument(
        "infile",
        nargs="?",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="Input text file (chrom\\tpos\\tref\\talt per line)",
    )
    parser.add_argument(
        "--pon-fasta", required=True, dest="pon_fasta",
        help="PON FASTA file (samtools faidx indexed)"
    )
    parser.add_argument(
        "--samtools-sif", required=True, dest="samtools_sif",
        help="Apptainer SIF for samtools"
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(name)s] %(levelname)s: %(message)s",
        stream=sys.stderr,
    )

    log.info("================================================================")
    log.info("START pon_mask_filter")
    log.info("input_file=%s", args.infile.name)
    log.info("pon_fasta=%s", args.pon_fasta)
    log.info("pon_fasta_size=%s bytes",
             os.path.getsize(args.pon_fasta) if os.path.exists(args.pon_fasta) else "N/A")
    log.info("samtools_sif=%s", args.samtools_sif)
    log.info("criterion: IUPAC match between alt allele and PON base → remove")
    log.info("================================================================")

    t0 = time.time()
    stats = filter_variants(args.infile, args.pon_fasta, args.samtools_sif)
    sys.stdout.flush()
    elapsed = time.time() - t0

    log.info("================================================================")
    log.info("RESULTS:")
    log.info("  variants_input=%d", stats["input"])
    log.info("  variants_kept=%d", stats["kept"])
    log.info("  variants_removed=%d (PON masked)", stats["removed"])
    log.info("  lines_skipped=%d (malformed)", stats["skipped"])
    log.info("  faidx_errors=%d", stats["faidx_errors"])
    if stats["input"] > 0:
        log.info("  pass_rate=%.1f%%", stats["kept"] / stats["input"] * 100)
    if stats.get("pon_code_distribution"):
        log.info("  PON base distribution:")
        for code, count in stats["pon_code_distribution"].items():
            log.info("    %s: %d", code, count)
    log.info("  elapsed=%.1f seconds", elapsed)
    log.info("END pon_mask_filter")
    log.info("================================================================")


if __name__ == "__main__":
    main()
