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
# Requirements: wget, bwa, samtools, gatk (or picard)
# RAM: BWA index requires ~32GB
# =============================================================================

set -uo pipefail

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
  echo "  Running samtools faidx..."
  samtools faidx "$REF_FINAL"
fi

if [[ -f "${REF_FINAL%.fasta}.dict" ]]; then
  echo "  [SKIP] .dict already exists"
else
  echo "  Creating sequence dictionary..."
  if command -v gatk &>/dev/null; then
    gatk CreateSequenceDictionary -R "$REF_FINAL"
  elif command -v picard &>/dev/null; then
    picard CreateSequenceDictionary R="$REF_FINAL" O="${REF_FINAL%.fasta}.dict"
  else
    samtools dict "$REF_FINAL" > "${REF_FINAL%.fasta}.dict"
  fi
fi

if [[ -f "${REF_FINAL}.sa" ]]; then
  echo "  [SKIP] BWA index already exists"
else
  echo "  Running bwa index (this takes ~1 hour, needs ~32GB RAM)..."
  bwa index "$REF_FINAL"
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

# ---- 4. File listing ------------------------------------------------------

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
