#!/usr/bin/env python3
"""Generate config/samples.tsv from a directory of paired FASTQ files.

Supports common naming conventions via auto-detection, custom regex, or
field-split mode. Handles multiple readgroups per sample automatically.

Usage examples:
    # Auto-detect naming convention
    python scripts/make_samples_tsv.py /data/fastq/

    # Custom regex with named groups
    python scripts/make_samples_tsv.py /data/fastq/ \\
        --pattern "(?P<sample_id>[^_]+)_(?P<rg>[^_]+)_R1"

    # Field-split: split by '_', use field[0] as sample_id, field[1] as rg
    python scripts/make_samples_tsv.py /data/fastq/ --split 0 1

    # Dry-run: preview without writing
    python scripts/make_samples_tsv.py /data/fastq/ --dry-run

    # Recursive directory scan
    python scripts/make_samples_tsv.py /data/fastq/ --recursive

    # Write to specific output file
    python scripts/make_samples_tsv.py /data/fastq/ -o config/samples.tsv
"""

import argparse
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# R1 / R2 marker definitions
# ---------------------------------------------------------------------------

# Each entry: (r1_pattern, r2_replacement)
# r1_pattern is matched against the full filename.
# r2_replacement is a substitution to derive the R2 filename.
R_MARKERS = [
    (r"_R1_001\.(fastq|fq)(\.gz)?$", "_R2_001"),  # Illumina BCL2FASTQ
    (r"_R1\.(fastq|fq)(\.gz)?$", "_R2"),  # simple _R1 / _R2
    (r"_1\.(fastq|fq)(\.gz)?$", "_2"),  # _1 / _2 (SRA, numeric)
    (r"\.R1\.(fastq|fq)(\.gz)?$", ".R2"),  # dot-separated .R1 / .R2
    (r"\.1\.(fastq|fq)(\.gz)?$", ".2"),  # dot-separated .1 / .2
]


# ---------------------------------------------------------------------------
# Built-in pattern registry
# ---------------------------------------------------------------------------
# Each pattern has:
#   name      — display name
#   regex     — must contain named groups 'sample_id' and optionally 'rg'
#               (if 'rg' group absent, rg is derived from sample_id)
#   rg_const  — if set, use this constant string as rg (overrides 'rg' group)

BUILTIN_PATTERNS = [
    {
        "name": "SRA",
        "desc": "SRR/ERR/DRR + numeric ID  (SRR123456_1.fastq.gz)",
        "regex": r"^(?P<sample_id>[SDE]RR\d+)_1\.(fastq|fq)(\.gz)?$",
        # rg = sample_id (each SRR accession is its own readgroup)
    },
    {
        "name": "Illumina-BCL2FASTQ",
        "desc": "Illumina standard  (SAMPLE_S1_L001_R1_001.fastq.gz)",
        "regex": r"^(?P<sample_id>.+?)_S\d+_(?P<rg>L\d{3})_R1_\d+\.(fastq|fq)(\.gz)?$",
    },
    {
        "name": "Lane-tagged",
        "desc": "Sample + lane tag  (SAMPLE_L001_R1.fastq.gz)",
        "regex": r"^(?P<sample_id>.+?)_(?P<rg>L\d{3})_R1\.(fastq|fq)(\.gz)?$",
    },
    {
        "name": "Run-tagged",
        "desc": "Sample + run/replicate  (SAMPLE_run1_R1.fastq.gz)",
        "regex": r"^(?P<sample_id>.+?)_(?P<rg>(?:run|rep|lib|RG)\w+)_R1\.(fastq|fq)(\.gz)?$",
    },
    {
        "name": "Simple-R1R2",
        "desc": "Simple paired  (SAMPLE_R1.fastq.gz)",
        "regex": r"^(?P<sample_id>.+?)_R1\.(fastq|fq)(\.gz)?$",
    },
    {
        "name": "Numeric",
        "desc": "Numeric suffix  (SAMPLE_1.fastq.gz) — non-SRA",
        "regex": r"^(?P<sample_id>.+?)_1\.(fastq|fq)(\.gz)?$",
    },
]


# ---------------------------------------------------------------------------
# Helper: find R1 files
# ---------------------------------------------------------------------------


def find_r1_files(directory: str, recursive: bool) -> list[Path]:
    """Return all FASTQ files that look like R1/read1."""
    root = Path(directory)
    if not root.is_dir():
        raise ValueError(f"Not a directory: {directory}")

    r1_files = []
    glob_fn = root.rglob if recursive else root.glob
    for p in sorted(glob_fn("*")):
        if not p.is_file():
            continue
        fname = p.name
        for r1_pat, _ in R_MARKERS:
            if re.search(r1_pat, fname):
                r1_files.append(p)
                break
    return r1_files


def find_r2(r1_path: Path) -> Path | None:
    """Derive and validate the R2 path from an R1 path."""
    fname = r1_path.name
    for r1_pat, r2_rep in R_MARKERS:
        m = re.search(r1_pat, fname)
        if m:
            # Reconstruct R2 filename by replacing the R1 marker
            ext = m.group(0)  # e.g. "_R1.fastq.gz"
            r2_ext = r2_rep + ext[ext.index(".") :]  # e.g. "_R2.fastq.gz"
            r2_name = fname[: m.start()] + r2_ext
            r2_path = r1_path.parent / r2_name
            return r2_path if r2_path.exists() else None
    return None


# ---------------------------------------------------------------------------
# Helper: extract sample_id and rg from filename
# ---------------------------------------------------------------------------


def _extract_via_regex(fname: str, regex: str) -> tuple[str, str | None] | None:
    """Return (sample_id, rg_or_None) if regex matches, else None."""
    m = re.match(regex, fname)
    if not m:
        return None
    gd = m.groupdict()
    sample_id = gd.get("sample_id")
    rg = gd.get("rg")  # may be None if group absent
    if not sample_id:
        return None
    return sample_id, rg


def extract_sample_rg(
    r1_path: Path,
    custom_pattern: str | None = None,
    split_fields: tuple[int, int] | None = None,
    split_delim: str = "_",
) -> tuple[str, str]:
    """Extract (sample_id, readgroup) from an R1 filename.

    Priority:
      1. custom_pattern (user-supplied regex)
      2. split_fields   (user-supplied field indices)
      3. built-in patterns (auto-detection)
      4. fallback: whole prefix = sample_id, rg = "RG1"
    """
    fname = r1_path.name

    # 1. Custom regex
    if custom_pattern:
        result = _extract_via_regex(fname, custom_pattern)
        if result:
            sample_id, rg = result
            return sample_id, rg if rg else sample_id
        # If custom pattern supplied but does not match, warn and fall through
        print(
            f"  WARNING: --pattern did not match '{fname}'; falling back to auto-detect.",
            file=sys.stderr,
        )

    # 2. Field split
    if split_fields is not None:
        # Strip R1 marker first to get clean base name
        base = fname
        for r1_pat, _ in R_MARKERS:
            base = re.sub(r1_pat, "", base)
        parts = base.split(split_delim)
        si_idx, rg_idx = split_fields
        try:
            sample_id = parts[si_idx]
            rg = parts[rg_idx]
            return sample_id, rg
        except IndexError:
            print(
                f"  WARNING: --split indices out of range for '{fname}'; "
                "falling back to auto-detect.",
                file=sys.stderr,
            )

    # 3. Built-in patterns
    for pat in BUILTIN_PATTERNS:
        result = _extract_via_regex(fname, pat["regex"])
        if result:
            sample_id, rg = result
            return sample_id, rg if rg else sample_id

    # 4. Fallback: strip R1 marker, use whole base as sample_id
    base = fname
    for r1_pat, _ in R_MARKERS:
        base = re.sub(r1_pat, "", base)
    return base, "RG1"


# ---------------------------------------------------------------------------
# Auto-detect dominant pattern across all R1 files
# ---------------------------------------------------------------------------


def detect_pattern(r1_files: list[Path]) -> str | None:
    """Return name of the built-in pattern that matches most R1 files."""
    counts: dict[str, int] = defaultdict(int)
    for r1 in r1_files:
        for pat in BUILTIN_PATTERNS:
            if re.match(pat["regex"], r1.name):
                counts[pat["name"]] += 1
                break
    if not counts:
        return None
    best = max(counts, key=lambda k: counts[k])
    coverage = counts[best] / len(r1_files)
    return best if coverage >= 0.5 else None


# ---------------------------------------------------------------------------
# Main logic
# ---------------------------------------------------------------------------


def build_table(
    directory: str,
    recursive: bool,
    custom_pattern: str | None,
    split_fields: tuple[int, int] | None,
    split_delim: str,
) -> list[dict[str, str]]:
    """Scan directory and return list of row dicts for samples.tsv."""
    r1_files = find_r1_files(directory, recursive)
    if not r1_files:
        raise RuntimeError(f"No R1 FASTQ files found in '{directory}'.")

    # Show detected pattern
    if custom_pattern:
        print(f"[make_samples_tsv] Using custom pattern: {custom_pattern}")
    elif split_fields:
        print(
            f"[make_samples_tsv] Using field-split: indices {split_fields}, delim='{split_delim}'"
        )
    else:
        detected = detect_pattern(r1_files)
        if detected:
            print(f"[make_samples_tsv] Auto-detected pattern: {detected}")
        else:
            print("[make_samples_tsv] No dominant pattern detected; using fallback (full prefix).")

    rows = []
    missing_r2 = []

    for r1 in r1_files:
        r2 = find_r2(r1)
        if r2 is None:
            missing_r2.append(r1)
            continue

        sample_id, rg = extract_sample_rg(
            r1,
            custom_pattern=custom_pattern,
            split_fields=split_fields,
            split_delim=split_delim,
        )
        rows.append(
            {
                "sample_id": sample_id,
                "readgroup": rg,
                "fq1": str(r1.resolve()),
                "fq2": str(r2.resolve()),
            }
        )

    if missing_r2:
        print(
            f"\n  WARNING: {len(missing_r2)} R1 file(s) have no matching R2 — skipped:",
            file=sys.stderr,
        )
        for p in missing_r2:
            print(f"    {p}", file=sys.stderr)

    return rows


def build_bam_row(alignment_path: Path, pattern: str | None = None) -> dict[str, str]:
    """Build a single bam-mode row from a pre-aligned ``.bam``/``.cram`` path.

    The ``bam`` value may be either a BAM or a CRAM path (the ingest rule
    normalizes both to CRAM downstream). ``sample_id`` is derived from the
    filename stem (extension stripped). If ``pattern`` is supplied, its named
    group ``sample_id`` is matched against the stem; on no match the full stem
    is used as a fallback.
    """
    stem = alignment_path.name
    for ext in (".bam", ".cram"):
        if stem.endswith(ext):
            stem = stem[: -len(ext)]
            break

    sample_id = stem
    if pattern:
        m = re.match(pattern, stem)
        if m and m.groupdict().get("sample_id"):
            sample_id = m.group("sample_id")
        else:
            print(
                f"  WARNING: --pattern did not match stem '{stem}'; using full stem as sample_id.",
                file=sys.stderr,
            )

    return {"sample_id": sample_id, "bam": str(alignment_path.resolve())}


def print_table(rows: list[dict[str, str]]) -> None:
    """Pretty-print the table to stdout."""
    header = f"{'sample_id':<20} {'readgroup':<20} {'fq1'}"
    print("\n" + "=" * 72)
    print(header)
    print("-" * 72)
    for r in rows:
        print(f"  {r['sample_id']:<18} {r['readgroup']:<18} {r['fq1']}")
        print(f"  {'':18} {'':18} {r['fq2']}")
    print("=" * 72)
    sample_count = len({r["sample_id"] for r in rows})
    print(f"  {len(rows)} readgroup(s) across {sample_count} sample(s)\n")


def _ordered_columns(rows: list[dict[str, str]]) -> list[str]:
    """Derive output columns from the row keys, sample_id always first.

    Schema-agnostic: works for both the fastq schema
    (sample_id, readgroup, fq1, fq2) and the bam schema (sample_id, bam).
    Column order follows first-row insertion order, with ``sample_id`` forced
    to the front so downstream readers can rely on it.
    """
    first = rows[0]
    cols = list(first.keys())
    if "sample_id" in cols:
        cols.remove("sample_id")
        cols.insert(0, "sample_id")
    return cols


def write_tsv(rows: list[dict[str, str]], output: str) -> None:
    """Write rows to a TSV file using the schema implied by the row keys."""
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    columns = _ordered_columns(rows)
    with open(output, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in columns) + "\n")
    print(f"[make_samples_tsv] Written: {output}  ({len(rows)} rows)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate samples.tsv from a directory of paired FASTQ files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "directory",
        nargs="?",
        help="Directory containing FASTQ files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="config/samples.tsv",
        help="Output TSV path (default: config/samples.tsv).",
    )
    parser.add_argument(
        "--pattern",
        default=None,
        metavar="REGEX",
        help=(
            "Custom regex with named groups (?P<sample_id>...) and "
            "optionally (?P<rg>...). Applied to R1 filenames."
        ),
    )
    parser.add_argument(
        "--split",
        nargs=2,
        type=int,
        metavar=("SAMPLE_IDX", "RG_IDX"),
        default=None,
        help=(
            "Field-split mode: split filename by --delim and use "
            "SAMPLE_IDX for sample_id, RG_IDX for readgroup. "
            "Example: --split 0 1"
        ),
    )
    parser.add_argument(
        "--delim",
        default="_",
        help="Delimiter for --split mode (default: '_').",
    )
    parser.add_argument(
        "--recursive",
        "-r",
        action="store_true",
        help="Scan subdirectories recursively.",
    )
    parser.add_argument(
        "--dry-run",
        "-n",
        action="store_true",
        dest="dry_run",
        help="Preview detected assignments without writing the output file.",
    )
    parser.add_argument(
        "--list-patterns",
        action="store_true",
        dest="list_patterns",
        help="Print built-in pattern registry and exit.",
    )
    args = parser.parse_args()

    if args.list_patterns:
        print("\nBuilt-in patterns (tried in order):")
        print("-" * 60)
        for i, pat in enumerate(BUILTIN_PATTERNS, 1):
            print(f"  {i}. {pat['name']}")
            print(f"     {pat['desc']}")
            print(f"     regex: {pat['regex']}")
        print()
        sys.exit(0)

    if not args.directory:
        parser.error("directory is required")

    split_fields = tuple(args.split) if args.split else None

    try:
        rows = build_table(
            directory=args.directory,
            recursive=args.recursive,
            custom_pattern=args.pattern,
            split_fields=split_fields,
            split_delim=args.delim,
        )
    except (ValueError, RuntimeError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    if not rows:
        print("ERROR: No valid FASTQ pairs found.", file=sys.stderr)
        sys.exit(1)

    print_table(rows)

    if args.dry_run:
        print("[make_samples_tsv] Dry-run mode: output file not written.")
    else:
        write_tsv(rows, args.output)


if __name__ == "__main__":
    main()
