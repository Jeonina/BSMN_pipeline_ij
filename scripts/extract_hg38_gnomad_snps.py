#!/usr/bin/env python3
"""Extract hg38 gnomAD SNPs above an AF threshold into a lookup TSV.

Generates the lookup table consumed by scripts/germline_filter.py from the
af-only gnomAD VCF (hg38 coordinates) that Mutect2 also uses as its
germline_resource.

Background (M-FIX-004): The legacy lookup `gnomAD.r2.1.1.AFover0.001.snps.txt.gz`
was hg19. Mutect2 calls hg38, so germline_filter cross-checked across builds
and let ~99% of common germline variants leak through. This script
regenerates the table directly from `af-only-gnomad.hg38.vcf.gz` so both
sides of the pipeline share the same coordinate system.

# @MX:ANCHOR: [AUTO] iter_passing_alleles is the canonical filter contract
# @MX:REASON: Output schema is consumed by scripts/germline_filter.load_germline_set;
# any change to row format breaks the downstream lookup key.

Output format (tab-separated, gzipped):
    chrom\tpos\tref\talt

Usage:
    python scripts/extract_hg38_gnomad_snps.py \\
        --input  resources/hg38/af-only-gnomad.hg38.vcf.gz \\
        --output resources/hg38/gnomAD.hg38.AFover0.001.snps.txt.gz \\
        --af-threshold 0.001
"""

from __future__ import annotations

import argparse
import gzip
import logging
import sys
import time
from collections.abc import Iterator
from pathlib import Path

log = logging.getLogger("extract_hg38_gnomad_snps")


def _parse_af_info(info: str) -> list[float] | None:
    """Extract the AF values from a VCF INFO column.

    Returns the list of parsed AF floats, or None when AF is absent or any
    value is malformed (e.g. `AF=.`, `AF=foo`). The whole record is treated
    as suspect on any parse failure — we never silently default to 0.
    """
    for field in info.split(";"):
        if not field.startswith("AF="):
            continue
        raw = field[3:]
        if not raw:
            return None
        out: list[float] = []
        for token in raw.split(","):
            try:
                out.append(float(token))
            except ValueError:
                return None
        return out
    return None


def iter_passing_alleles(
    line: str, af_threshold: float, strip_chr: bool
) -> Iterator[tuple[str, str, str, str]]:
    """Yield (chrom, pos, ref, alt) tuples for SNV alts passing the threshold.

    Handles multiallelic records by splitting ALT on `,` and matching each
    alt to its index-aligned AF value. Skips:
      * indels (ref length != 1, or alt length != 1)
      * AF <= af_threshold
      * records with malformed/missing AF
      * records where len(alts) != len(AF)
    """
    line = line.rstrip("\n")
    if not line or line.startswith("#"):
        return
    parts = line.split("\t")
    if len(parts) < 8:
        return
    chrom, pos, _id, ref, alt_field, _qual, _filter, info = parts[:8]

    if len(ref) != 1:
        return

    afs = _parse_af_info(info)
    if afs is None:
        return

    alts = alt_field.split(",")
    if len(alts) != len(afs):
        return

    emit_chrom = chrom[3:] if strip_chr and chrom.startswith("chr") else chrom

    for alt, af in zip(alts, afs, strict=True):
        if len(alt) != 1:
            continue
        if af <= af_threshold:
            continue
        yield emit_chrom, pos, ref, alt


def extract(
    input_path: Path,
    output_path: Path,
    af_threshold: float,
    strip_chr: bool,
) -> dict[str, int]:
    """Stream `input_path` (gzipped VCF) into `output_path` (gzipped TSV).

    Returns counters for logging: total lines read, biallelic SNV records
    that emitted output, multiallelic SNV records, alleles passing the
    threshold (== output rows), and skipped records.
    """
    if not input_path.exists():
        raise FileNotFoundError(f"input VCF not found: {input_path}")

    stats: dict[str, int] = {
        "lines": 0,
        "records": 0,
        "biallelic_emitted": 0,
        "multiallelic_emitted": 0,
        "alleles_emitted": 0,
        "skipped": 0,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(input_path, "rt") as src, gzip.open(
        output_path, "wt", compresslevel=6
    ) as dst:
        for line in src:
            stats["lines"] += 1
            if not line or line.startswith("#"):
                continue
            stats["records"] += 1
            rows = list(iter_passing_alleles(line, af_threshold, strip_chr))
            if not rows:
                stats["skipped"] += 1
                continue
            if len(rows) == 1:
                stats["biallelic_emitted"] += 1
            else:
                stats["multiallelic_emitted"] += 1
            stats["alleles_emitted"] += len(rows)
            for chrom, pos, ref, alt in rows:
                dst.write(f"{chrom}\t{pos}\t{ref}\t{alt}\n")

    return stats


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Extract hg38 gnomAD SNPs above an AF threshold into the lookup "
            "table consumed by scripts/germline_filter.py."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Path to gzipped af-only-gnomad.hg38.vcf.gz",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Path to gzipped output TSV (chrom\\tpos\\tref\\talt)",
    )
    parser.add_argument(
        "--af-threshold",
        type=float,
        default=0.001,
        help="Minimum (exclusive) allele frequency to retain (default: 0.001)",
    )
    parser.add_argument(
        "--no-strip-chr",
        action="store_true",
        help=(
            "Keep the 'chr' prefix in the output chrom column. Default behaviour "
            "strips it to match the legacy gnomAD lookup format."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )
    args = _build_parser().parse_args(argv)

    t0 = time.perf_counter()
    try:
        stats = extract(
            input_path=args.input,
            output_path=args.output,
            af_threshold=args.af_threshold,
            strip_chr=not args.no_strip_chr,
        )
    except FileNotFoundError as exc:
        log.error("%s", exc)
        return 2
    except OSError as exc:
        log.error("I/O error reading/writing VCF: %s", exc)
        return 3

    elapsed = time.perf_counter() - t0
    log.info(
        "lines=%d records=%d biallelic_emitted=%d multiallelic_emitted=%d "
        "alleles_emitted=%d skipped=%d elapsed=%.2fs",
        stats["lines"],
        stats["records"],
        stats["biallelic_emitted"],
        stats["multiallelic_emitted"],
        stats["alleles_emitted"],
        stats["skipped"],
        elapsed,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
