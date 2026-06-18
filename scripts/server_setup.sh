#!/usr/bin/env bash
# =============================================================================
# server_setup.sh — BSMN pipeline server setup
#
# Run this script ON the server after cloning the repo.
#
# Usage:
#   bash scripts/server_setup.sh
#
# Requirements on server:
#   - git, apptainer, conda (or mamba), wget, samtools, lftp
#   - Internet access (for resource download and container pull)
#   - ~60 GB free disk space for resources/hg38/
#
# Resources are obtained in two ways:
#   1. Public files — downloaded automatically by download_and_index_hg38.sh
#   2. Custom files (gnomAD, PON) — assembled from downloads/ (git-tracked split parts)
# =============================================================================
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"
echo "[setup] Working directory: $REPO_DIR"

# ---------------------------------------------------------------------------
# 1. Conda environment
# ---------------------------------------------------------------------------
echo ""
echo "[1/5] Setting up conda environment..."

if conda env list | grep -q "^bp "; then
    echo "  → conda env 'bp' already exists, skipping create"
else
    conda env create -f environment.yml
    echo "  → conda env 'bp' created"
fi

# Activate for subsequent commands
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bp

# ---------------------------------------------------------------------------
# 2. Pull Apptainer containers
# ---------------------------------------------------------------------------
echo ""
echo "[2/5] Pulling Apptainer containers..."
echo "  (This may take 10–30 minutes depending on network speed)"

python3 scripts/prepare_containers.py config/containers.yaml
echo "  → All containers ready"

# ---------------------------------------------------------------------------
# 3. Download public reference resources
# ---------------------------------------------------------------------------
echo ""
echo "[3/5] Downloading public reference resources..."
echo "  (This may take 1–3 hours depending on network speed)"
echo "  Requires: wget, samtools, lftp"

mkdir -p resources/hg38
bash scripts/download_and_index_hg38.sh resources/hg38
echo "  → Public resources downloaded"

# ---------------------------------------------------------------------------
# 4. Assemble split resource files + clone MosaicForecast model
# ---------------------------------------------------------------------------
echo ""
echo "[4/5] Assembling custom resources..."

# NOTE: the gnomAD germline lookup (gnomAD.hg38.AFover0.001.snps.txt.gz) is built
# by download_and_index_hg38.sh Step 6 from the hg38 af-only VCF — do NOT assemble
# the legacy hg19 gnomAD.r2.1.1 lookup here (M-FIX-004: hg19 coords let ~99% of
# common germline variants leak through germline_filter).

# PON: cat parts → gunzip → resources/hg38/
PON_OUT="resources/hg38/PON.q20q20.05.5.fa"
PON_PARTS=(
    downloads/PON.q20q20.05.5.fa.gz.partaa
    downloads/PON.q20q20.05.5.fa.gz.partab
    downloads/PON.q20q20.05.5.fa.gz.partac
    downloads/PON.q20q20.05.5.fa.gz.partad
)
if [[ ! -f "$PON_OUT" ]]; then
    echo "  → Assembling PON..."
    cat "${PON_PARTS[@]}" | gunzip -c > "$PON_OUT"
    samtools faidx "$PON_OUT"
    echo "  ✓ $PON_OUT"
else
    echo "  ✓ $PON_OUT (already assembled)"
fi

# MosaicForecast: clone repo for the trained RF model (default E-step model,
# config filtering.mosaicforecast.model). The container ships the MF scripts +
# k24 bigwig but NOT the models_trained/*.rds files.
MF_DIR="resources/MosaicForecast"
MF_MODEL="$MF_DIR/models_trained/250xRFmodel_addRMSK_Refine.rds"
if [[ -f "$MF_MODEL" ]]; then
    echo "  ✓ $MF_MODEL (already present)"
else
    echo "  → Cloning MosaicForecast (for trained RF model)..."
    rm -rf "$MF_DIR"
    git clone --depth 1 https://github.com/parklab/MosaicForecast.git "$MF_DIR"
    if [[ -f "$MF_MODEL" ]]; then
        echo "  ✓ $MF_MODEL"
    else
        echo "  ✗ MosaicForecast clone did not contain $MF_MODEL"
    fi
fi

# ---------------------------------------------------------------------------
# 5. Verify resources
# ---------------------------------------------------------------------------
echo ""
echo "[5/5] Verifying resources..."

REQUIRED_FILES=(
    "resources/hg38/Homo_sapiens_assembly38.fasta"
    "resources/hg38/Homo_sapiens_assembly38.fasta.bwt"
    "resources/hg38/Homo_sapiens_assembly38.fasta.fai"
    "resources/hg38/Homo_sapiens_assembly38.dict"
    "resources/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
    "resources/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    "resources/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    "resources/hg38/af-only-gnomad.hg38.vcf.gz"
    "resources/hg38/small_exac_common_3.hg38.vcf.gz"
    "resources/hg38/1KG.20160622.strict_mask.hg38_GRCh38.bed"
    "resources/hg38/gnomAD.hg38.AFover0.001.snps.txt.gz"
    "resources/hg38/PON.q20q20.05.5.fa"
    "resources/hg38/PON.q20q20.05.5.fa.fai"
    "resources/MosaicForecast/models_trained/250xRFmodel_addRMSK_Refine.rds"
)

ALL_OK=true
for f in "${REQUIRED_FILES[@]}"; do
    if [[ -f "$f" ]]; then
        echo "  ✓ $f"
    else
        echo "  ✗ MISSING: $f"
        ALL_OK=false
    fi
done

if [[ "$ALL_OK" == false ]]; then
    echo ""
    echo "[ERROR] Some resource files are missing."
    echo "  Check that Step 3 (download) and Step 4 (assemble) completed without errors."
    echo "  For custom files (gnomAD, PON), ensure downloads/ split parts are present:"
    echo "    ls downloads/"
    echo "  Then re-run this script."
    exit 1
fi

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
echo ""
echo "=========================================="
echo " Setup complete. Run the pipeline:"
echo ""
echo "  conda activate bp"
echo "  cd $REPO_DIR"
echo ""
echo "  # Mapping only:"
echo "  snakemake --snakefile workflow/Snakefile --cores 16"
echo ""
echo "  # Through calling:"
echo "  snakemake --snakefile workflow/Snakefile --config stage=calling --cores 16"
echo ""
echo "  # Full pipeline (mapping + calling + filtering):"
echo "  snakemake --snakefile workflow/Snakefile --config stage=filtering --cores 16"
echo "=========================================="
