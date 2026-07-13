#!/usr/bin/env bash
# =============================================================================
# Download hg38 resources & build BWA index
# Based on BSMN pipeline (Broad GATK bundle)
# Sources: NCBI, EBI (Google Cloud Storage blocked)
#
# Usage:
#   bash download_and_index_hg38.sh [output_dir]
#   (default output_dir: ./resources_hg38)
#
# Requirements: wget, lftp, apptainer, python3 (with pyyaml)
# Bio tools (bwa/samtools/gatk) run inside Apptainer containers (config/containers.yaml).
# RAM: BWA index requires ~32GB
# =============================================================================

set -uo pipefail

# Resolve the scripts/ dir BEFORE cd into OUTDIR (BASH_SOURCE is relative to the
# caller's CWD; computing it after the cd below would break — see Step 6).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Resolve Apptainer SIF paths from config/containers.yaml (single source of truth).
# Needs pyyaml (provided by the project .venv / requirements.txt).
_sif() {
  python3 -c "import yaml,sys; print(yaml.safe_load(open(sys.argv[1]))[sys.argv[2]]['sif'])" \
    "$REPO_ROOT/config/containers.yaml" "$1"
}
BWA_SIF="$REPO_ROOT/$(_sif bwa)"
SAMTOOLS_SIF="$REPO_ROOT/$(_sif samtools)"
GATK_SIF="$REPO_ROOT/$(_sif gatk)"

OUTDIR=${1:-./resources_hg38}
mkdir -p "$OUTDIR"
cd "$OUTDIR"

SECONDS=0

echo "============================================="
echo " Output directory: $(pwd)"
echo "============================================="

# ---- 1. Reference FASTA ---------------------------------------------------

echo ""
echo "[Step 1] Downloading reference FASTA..."

NCBI_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids"
REF_GZ="GCA_000001405.15_GRCh38_full_analysis_set.fna.gz"
REF_FINAL="Homo_sapiens_assembly38.fasta"

if [[ -f "$REF_FINAL" ]]; then
  echo "  [SKIP] $REF_FINAL already exists"
else
  echo "  Downloading $REF_GZ ..."
  wget -c -q --show-progress "$NCBI_BASE/$REF_GZ"
  echo "  Decompressing..."
  gunzip -c "$REF_GZ" > "$REF_FINAL"
  rm -f "$REF_GZ"
fi

echo "[Step 1] Done."

# ---- 2. Build index files -------------------------------------------------

echo ""
echo "[Step 2] Building index files..."

if [[ -f "${REF_FINAL}.fai" ]]; then
  echo "  [SKIP] .fai already exists"
else
  echo "  Running samtools faidx (via Apptainer)..."
  apptainer exec "$SAMTOOLS_SIF" samtools faidx "$REF_FINAL"
fi

if [[ -f "${REF_FINAL%.fasta}.dict" ]]; then
  echo "  [SKIP] .dict already exists"
else
  echo "  Creating sequence dictionary (via Apptainer)..."
  apptainer exec "$GATK_SIF" gatk CreateSequenceDictionary -R "$REF_FINAL"
fi

if [[ -f "${REF_FINAL}.sa" ]]; then
  echo "  [SKIP] BWA index already exists"
else
  echo "  Running bwa index (via Apptainer; this takes ~1 hour, needs ~32GB RAM)..."
  apptainer exec "$BWA_SIF" bwa index "$REF_FINAL"
fi

echo "[Step 2] Done."

# ---- 3. Known sites VCFs (for BQSR) ---------------------------------------

echo ""
echo "[Step 3] Downloading known sites VCFs..."

FAIL_COUNT=0

download_vcf() {
  local dest="$1"
  local url="$2"
  if [[ -f "$dest" ]]; then
    echo "  [SKIP] $dest already exists"
  else
    echo "  Downloading $dest ..."
    if wget -c -q --show-progress -O "$dest" "$url"; then
      echo "  [OK] $dest"
    else
      echo "  [FAIL] $dest — URL: $url"
      rm -f "$dest"
      FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
  fi
}

# 3a. dbSNP - NCBI
DBSNP_URL="https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK"
download_vcf "Homo_sapiens_assembly38.dbsnp138.vcf.gz"     "$DBSNP_URL/00-All.vcf.gz"
download_vcf "Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi" "$DBSNP_URL/00-All.vcf.gz.tbi"

# 3b. Mills indels - EBI (try multiple paths)
MILLS="Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
EBI_1="https://ftp.ebi.ac.uk/pub/databases/1000genomes/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources"
EBI_2="https://ftp.ebi.ac.uk/pub/databases/1000genomes/ftp/technical/reference/GRCh38_reference_genome"
for suffix in "" ".tbi"; do
  if [[ -f "${MILLS}${suffix}" ]]; then
    echo "  [SKIP] ${MILLS}${suffix} already exists"
  else
    echo "  Downloading ${MILLS}${suffix} (trying EBI path 1)..."
    if wget -c -q --show-progress -O "${MILLS}${suffix}" "$EBI_1/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz${suffix}" 2>/dev/null; then
      echo "  [OK] ${MILLS}${suffix}"
    else
      echo "  EBI path 1 failed, trying path 2..."
      if wget -c -q --show-progress -O "${MILLS}${suffix}" "$EBI_2/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz${suffix}" 2>/dev/null; then
        echo "  [OK] ${MILLS}${suffix}"
      else
        echo "  [FAIL] ${MILLS}${suffix} — could not download from EBI"
        rm -f "${MILLS}${suffix}"
        FAIL_COUNT=$((FAIL_COUNT + 1))
      fi
    fi
  fi
done

# 3c. 1000G SNPs - EBI
KG="1000G_phase1.snps.high_confidence.hg38.vcf.gz"
for suffix in "" ".tbi"; do
  if [[ -f "${KG}${suffix}" ]]; then
    echo "  [SKIP] ${KG}${suffix} already exists"
  else
    echo "  Downloading ${KG}${suffix} (trying EBI)..."
    if wget -c -q --show-progress -O "${KG}${suffix}" "$EBI_1/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapped.20150424.genotypes.vcf.gz${suffix}" 2>/dev/null; then
      echo "  [OK] ${KG}${suffix}"
    else
      echo "  [FAIL] ${KG}${suffix} — could not download from EBI"
      rm -f "${KG}${suffix}"
      FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
  fi
done

if [[ "$FAIL_COUNT" -gt 0 ]]; then
  echo ""
  echo "  [WARNING] $FAIL_COUNT file(s) failed to download. See above for details."
else
  echo ""
fi
echo "[Step 3] Done."

# ---- 4. Contamination resource (small_exac_common) -----------------------

echo ""
echo "[Step 4] Downloading contamination resource (small_exac_common)..."

EXAC="small_exac_common_3.hg38.vcf.gz"
BROAD_FTP="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2"

if [[ -f "$EXAC" && -f "${EXAC}.tbi" ]]; then
  echo "  [SKIP] $EXAC already exists"
else
  echo "  Downloading $EXAC via FTP (requires lftp)..."
  lftp -c "open ${BROAD_FTP}; get ${EXAC} -o ${EXAC}; get ${EXAC}.tbi -o ${EXAC}.tbi"
fi

echo "[Step 4] Done."

# ---- 5. 1KG strict mask (mappability filter) ------------------------------

echo ""
echo "[Step 5] Downloading 1KG strict mask (BED)..."

MASK_OUT="1KG.20160622.strict_mask.hg38_GRCh38.bed"
EBI_MASK="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/StrictMask/20160622.allChr.mask.bed"

if [[ -f "$MASK_OUT" ]]; then
  echo "  [SKIP] $MASK_OUT already exists"
else
  echo "  Downloading 1KG strict mask from EBI..."
  wget -c -q --show-progress -O "$MASK_OUT" "$EBI_MASK"
fi

echo "[Step 5] Done."

# ---- 6. gnomAD hg38 SNP lookup (M-FIX-004) --------------------------------
#
# The germline_filter rule consumes a flat lookup table
# (chrom\tpos\tref\talt, gzipped) derived from af-only-gnomad.hg38.vcf.gz —
# the same hg38 resource Mutect2 uses as --germline-resource. Generating the
# lookup from the hg38 VCF guarantees the lookup coordinates match Mutect2's
# calls. Prior to M-FIX-004 the pipeline shipped an hg19-coordinate lookup
# (gnomAD.r2.1.1.AFover0.001.snps.txt.gz), which let ~99% of common germline
# variants leak through germline_filter.

echo ""
echo "[Step 6] gnomAD hg38: download af-only VCF + build SNP lookup..."

GNOMAD_VCF="af-only-gnomad.hg38.vcf.gz"
GNOMAD_LOOKUP="gnomAD.hg38.AFover0.001.snps.txt.gz"
# SCRIPT_DIR is captured at the top of this script (before cd "$OUTDIR").
BROAD_FTP="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/Mutect2"

# 6a. af-only-gnomad.hg38.vcf.gz — this is ALSO Mutect2's --germline-resource
#     (config calling.germline_resource). Building the germline lookup from this
#     same hg38 VCF guarantees the lookup coordinates match Mutect2's calls.
if [[ -f "$GNOMAD_VCF" && -f "${GNOMAD_VCF}.tbi" ]]; then
  echo "  [SKIP] $GNOMAD_VCF already exists"
else
  echo "  Downloading $GNOMAD_VCF via FTP (requires lftp)..."
  lftp -c "open ${BROAD_FTP}; get ${GNOMAD_VCF} -o ${GNOMAD_VCF}; get ${GNOMAD_VCF}.tbi -o ${GNOMAD_VCF}.tbi"
fi

# 6b. Flat SNP lookup for germline_filter (M-FIX-004: MUST be hg38 coords)
if [[ ! -f "$GNOMAD_VCF" ]]; then
  echo "  [FAIL] $GNOMAD_VCF missing after download — cannot build lookup"
  FAIL_COUNT=$((FAIL_COUNT + 1))
elif [[ -f "$GNOMAD_LOOKUP" ]]; then
  echo "  [SKIP] $GNOMAD_LOOKUP already exists"
else
  echo "  Building $GNOMAD_LOOKUP from $GNOMAD_VCF ..."
  python3 "${SCRIPT_DIR}/extract_hg38_gnomad_snps.py" \
    --input "$GNOMAD_VCF" \
    --output "$GNOMAD_LOOKUP" \
    --af-threshold 0.001
fi

echo "[Step 6] Done."

# ---- 7. File listing ------------------------------------------------------

echo ""
echo "============================================="
echo " Download complete. File listing:"
echo "============================================="
ls -lh

echo ""
echo "Generating md5 checksums..."
md5sum *.fasta *.gz > md5sums.txt 2>/dev/null || md5 *.fasta *.gz > md5sums.txt 2>/dev/null
echo "Saved to md5sums.txt"

elapsed=$SECONDS
printf "\n>> Total %d hours, %d minutes, %d seconds elapsed.\n" \
  $((elapsed / 3600)) $((elapsed % 3600 / 60)) $((elapsed % 60))

echo ""
echo "============================================="
echo " Next: copy this directory to your target server"
echo "   scp -r $(pwd)/ user@server:~/BSMN_pipeline_ij/resources/hg38/"
echo "============================================="
