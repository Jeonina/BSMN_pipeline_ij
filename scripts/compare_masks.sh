#!/usr/bin/env bash
# =============================================================================
# compare_masks.sh — quantify the overlap between the 1KG accessibility mask
# and the GIAB HG002 mosaic confident region.
#
# Answers the question: are the two region sets independent, overlapping, or is
# one an extension (superset) of the other?
#
#   - 1KG strict mask (accessibility) : population mappability mask, used as the
#     production accessibility filter (Step B). BED, 0-based half-open.
#   - GIAB HG002 mosaic confident BED : benchmark-scoring region for HG002
#     (a 2.45 Gbp subset of the GIAB germline benchmark regions, autosomes only).
#
# The two are built by different consortia from different evidence, so any
# overlap is incidental, not by construction. This script measures it.
#
# Requirements (on the data server): bedtools, awk, sort, curl/wget, gzip.
#
# Usage:
#   scripts/compare_masks.sh \
#       [KG_BED] [GIAB_BED] [OUTDIR]
#
# Defaults:
#   KG_BED   = resources/hg38/1KG.20160622.strict_mask.hg38_GRCh38.bed
#   GIAB_BED = resources/hg38/HG002_GRCh38_MosaicSNVv1.1.bed
#              (auto-downloaded from the GIAB FTP if missing)
#   OUTDIR   = reports/mask_overlap
#
# Output: a summary table (bp, %, Jaccard) printed to stdout and written to
#   OUTDIR/mask_overlap_summary.txt, plus the intersection/only BED files.
# =============================================================================
set -euo pipefail

KG_BED="${1:-resources/hg38/1KG.20160622.strict_mask.hg38_GRCh38.bed}"
GIAB_BED="${2:-resources/hg38/HG002_GRCh38_MosaicSNVv1.1.bed}"
OUTDIR="${3:-reports/mask_overlap}"

GIAB_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/mosaic_v1.10/GRCh38/SNV/HG002_GRCh38_MosaicSNVv1.1.bed"

command -v bedtools >/dev/null || { echo "ERROR: bedtools not found in PATH" >&2; exit 1; }
mkdir -p "$OUTDIR"

# --- Fetch GIAB confident BED if absent ------------------------------------
if [[ ! -s "$GIAB_BED" ]]; then
    echo "[info] GIAB confident BED not found at $GIAB_BED — downloading..." >&2
    mkdir -p "$(dirname "$GIAB_BED")"
    if command -v curl >/dev/null; then
        curl -fSL "$GIAB_URL" -o "$GIAB_BED"
    else
        wget -O "$GIAB_BED" "$GIAB_URL"
    fi
fi
[[ -s "$KG_BED" ]] || { echo "ERROR: 1KG mask BED not found at $KG_BED" >&2; exit 1; }

# --- Normalize: force chr-prefix, keep autosomes only, strip extra cols ------
# GIAB mosaic confident is autosomes only, so we compare on chr1..chr22 to be
# fair. Both inputs are coerced to chr-prefixed naming.
norm() {
    # $1 = input bed (may be .gz)
    local f="$1"
    local reader="cat"
    [[ "$f" == *.gz ]] && reader="zcat"
    $reader "$f" \
        | grep -v -E '^(#|track|browser)' \
        | awk 'BEGIN{OFS="\t"}
               {
                 c=$1;
                 if (c !~ /^chr/) c="chr"c;     # add chr prefix if missing
                 if (c ~ /^chr([1-9]|1[0-9]|2[0-2])$/) print c,$2,$3
               }' \
        | sort -k1,1 -k2,2n
}

KG_N="$OUTDIR/kg.autosomes.bed"
GIAB_N="$OUTDIR/giab.autosomes.bed"
norm "$KG_BED"   | bedtools merge -i - > "$KG_N"
norm "$GIAB_BED" | bedtools merge -i - > "$GIAB_N"

bp() { awk '{s+=$3-$2} END{printf "%.0f", s+0}' "$1"; }

A=$(bp "$KG_N")            # 1KG accessible bp (autosomes)
B=$(bp "$GIAB_N")          # GIAB confident bp (autosomes)

bedtools intersect -a "$KG_N" -b "$GIAB_N" > "$OUTDIR/intersection.bed"
bedtools subtract  -a "$KG_N" -b "$GIAB_N" > "$OUTDIR/kg_only.bed"
bedtools subtract  -a "$GIAB_N" -b "$KG_N" > "$OUTDIR/giab_only.bed"

I=$(bp "$OUTDIR/intersection.bed")     # intersection bp
KG_ONLY=$(bp "$OUTDIR/kg_only.bed")    # in 1KG, not GIAB
GIAB_ONLY=$(bp "$OUTDIR/giab_only.bed")# in GIAB, not 1KG

# Jaccard (bedtools needs identically-sorted inputs — both already sorted+merged)
JACC_LINE=$(bedtools jaccard -a "$KG_N" -b "$GIAB_N" | tail -n1)

SUMMARY="$OUTDIR/mask_overlap_summary.txt"
{
echo "=============================================================="
echo " 1KG accessibility  vs  GIAB HG002 mosaic confident  (autosomes)"
echo "=============================================================="
echo "1KG mask BED        : $KG_BED"
echo "GIAB confident BED  : $GIAB_BED"
echo "--------------------------------------------------------------"
awk -v A="$A" -v B="$B" -v I="$I" -v KO="$KG_ONLY" -v GO="$GIAB_ONLY" 'BEGIN{
  printf "1KG accessible (A)        : %15d bp\n", A
  printf "GIAB confident (B)        : %15d bp\n", B
  printf "Intersection (A∩B)        : %15d bp\n", I
  printf "1KG-only   (A∖B)          : %15d bp\n", KO
  printf "GIAB-only  (B∖A)          : %15d bp\n", GO
  printf "--------------------------------------------------------------\n"
  if (A>0) printf "Fraction of 1KG  in GIAB  : %6.2f %%   (A∩B / A)\n", 100*I/A
  if (B>0) printf "Fraction of GIAB in 1KG   : %6.2f %%   (A∩B / B)\n", 100*I/B
}'
echo "--------------------------------------------------------------"
echo "bedtools jaccard (intersection union jaccard n_intersections):"
echo "  $JACC_LINE"
echo "=============================================================="
echo "Interpretation:"
echo "  * Both containment fractions high (e.g. >90%) AND neither =100%"
echo "    => heavily OVERLAPPING but INDEPENDENT (neither is a subset)."
echo "  * One fraction ~100% => that set is (nearly) a SUBSET of the other"
echo "    (i.e. the other is an extension/superset)."
echo "  * Both fractions low => the two are largely DISJOINT."
echo "=============================================================="
} | tee "$SUMMARY"

echo "" >&2
echo "[done] summary -> $SUMMARY" >&2
echo "[done] BEDs    -> $OUTDIR/{intersection,kg_only,giab_only}.bed" >&2
