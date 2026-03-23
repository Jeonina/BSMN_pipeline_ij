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
import subprocess
import sys
from typing import Set

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
    except subprocess.CalledProcessError:
        return "?"


def filter_variants(txt_file, pon_fasta: str, samtools_sif: str, outfile=None) -> None:
    """Filter text-format variant file by PON IUPAC mask.

    Input format: chrom  pos  ref  alt  (tab-separated)
    Keeps variants where pon_passes(alt, PON_base) is True.
    """
    if outfile is None:
        outfile = sys.stdout

    for line in txt_file:
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 4:
            continue
        chrom, pos, _ref, alt = parts[0], parts[1], parts[2], parts[3]
        pon_base = query_fasta(pon_fasta, chrom, pos, samtools_sif)
        if pon_passes(alt, pon_base):
            outfile.write(line)


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

    filter_variants(args.infile, args.pon_fasta, args.samtools_sif)
    sys.stdout.flush()


if __name__ == "__main__":
    main()
