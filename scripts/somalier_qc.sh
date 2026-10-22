#!/usr/bin/env bash
# Sample QC before burden analysis with somalier: sex, identity/relatedness, genetic ancestry.
#
# Germline-WES habit (sex / kinship / population PCA) carried over to the somatic cohort:
#   relate   -> sex (X het, Y depth) vs metadata, and pairwise relatedness. Here the two
#               regions of one donor MUST be ~1 (same person); ~1 between different donors
#               = swap/duplicate. Kinship between relatives does not inflate somatic
#               burden (post-zygotic), so relatives are only noted, not removed.
#   ancestry -> 1000 Genomes super-population per sample (self-reported Race in the
#               Phenotype sheet is checked against it; Dx groups differ in Race).
#
# Usage (VM100, pipeline directory; the CRAM tree under results/):
#   Counting uses the bcftools container (APPTAINER_BIND defaults to /mnt/data,/tmp).
#   Re-running skips finished extracts.
#   bash ~/somalier_qc.sh ~/sample_info.tsv /mnt/data/bsmn/somalier 24 6   # jobs, autosomal stride
set -euo pipefail

INFO=${1:?sample_info.tsv (sample_id, br_region, sex)}
OUT=${2:?output dir}
JOBS=${3:-8}
STRIDE=${4:-6}     # count every 6th autosomal site (~2.9k): I/O is ~3 MB of CRAM per site
REF=${REF:-resources/hg38/Homo_sapiens_assembly38.fasta}
VER=v0.3.5
mkdir -p "$OUT/extracted" "$OUT/bin"
cd_pipe=$(pwd)

# --- one-time downloads (versions pinned) -----------------------------------------
S=$OUT/bin/somalier
[ -x "$S" ] || { curl -fsSL -o "$S" https://github.com/brentp/somalier/releases/download/$VER/somalier; chmod +x "$S"; }
[ -s "$OUT/sites.hg38.vcf.gz" ] || curl -fsSL -o "$OUT/sites.hg38.vcf.gz" \
    https://github.com/brentp/somalier/files/3412456/sites.hg38.vcf.gz
[ -s "$OUT/ancestry-labels-1kg.tsv" ] || curl -fsSL -o "$OUT/ancestry-labels-1kg.tsv" \
    https://raw.githubusercontent.com/brentp/somalier/master/scripts/ancestry-labels-1kg.tsv
[ -d "$OUT/1kg-somalier" ] || { curl -fsSL -o "$OUT/1kg.somalier.tar.gz" \
    "https://zenodo.org/record/3479773/files/1kg.somalier.tar.gz?download=1"
    tar -xzf "$OUT/1kg.somalier.tar.gz" -C "$OUT"; }
echo "somalier $VER"   # `somalier --version` exits non-zero, which set -e would treat as failure

# --- extract: one .somalier per CRAM (sample name = D-number) ----------------------
# `somalier extract` took ~5 h per ~100x CRAM over NFS: each site costs ~3 MB of CRAM
# reads and the storage tops out near 70 MB/s for any tool. So bcftools counts only a
# stride of the autosomal sites (plus all X/Y) and somalier_counts.py writes the
# .somalier binary (uncounted sites = 0/0/0 = no data for somalier).
export APPTAINER_BIND=${APPTAINER_BIND:-/mnt/data,/tmp}
python "$(dirname "$0")/somalier_counts.py" --sample-info "$INFO" --results results \
    --sites "$OUT/sites.hg38.vcf.gz" --ref "$REF" \
    --bcftools "apptainer exec containers/bcftools_1.17.sif bcftools" \
    --out-dir "$OUT/extracted" --jobs "$JOBS" --auto-stride "$STRIDE" || echo "WARNING: some samples failed; relating the rest"
echo "extracted: $(ls "$OUT"/extracted/*.somalier | wc -l)"

# --- relate: ped carries the metadata sex (1=M, 2=F); family = donor ----------------
awk -F'\t' 'NR>1{split($2,a,"_"); print a[1]"\t"$1"\t0\t0\t"($3=="M"?1:2)"\t-9"}' "$INFO" > "$OUT/cohort.ped"
"$S" relate --ped "$OUT/cohort.ped" -o "$OUT/relate" "$OUT"/extracted/*.somalier

# --- ancestry against 1000 Genomes ---------------------------------------------------
"$S" ancestry --labels "$OUT/ancestry-labels-1kg.tsv" -o "$OUT/ancestry" \
    "$OUT"/1kg-somalier/*.somalier ++ "$OUT"/extracted/*.somalier
ls -1 "$OUT"/relate.* "$OUT"/ancestry.*
