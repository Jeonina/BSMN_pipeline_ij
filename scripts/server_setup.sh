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
#   - git, apptainer, conda (or mamba)
#   - Internet access (for container pull)
#   - resources/ already transferred via rsync (see sync_resources.sh)
# =============================================================================
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"
echo "[setup] Working directory: $REPO_DIR"

# ---------------------------------------------------------------------------
# 1. Conda environment
# ---------------------------------------------------------------------------
echo ""
echo "[1/3] Setting up conda environment..."

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
echo "[2/3] Pulling Apptainer containers..."
echo "  (This may take 10–30 minutes depending on network speed)"

python3 scripts/prepare_containers.py config/containers.yaml
echo "  → All containers ready"

# ---------------------------------------------------------------------------
# 3. Assemble split resource files
# ---------------------------------------------------------------------------
echo ""
echo "[3/4] Assembling split resource files..."

# gnomAD: cat parts → resources/hg38/
GNOMAD_OUT="resources/hg38/gnomAD.r2.1.1.AFover0.001.snps.txt.gz"
GNOMAD_PARTS=(
    downloads/gnomAD.r2.1.1.AFover0.001.snps.txt.gz.partaa
    downloads/gnomAD.r2.1.1.AFover0.001.snps.txt.gz.partab
    downloads/gnomAD.r2.1.1.AFover0.001.snps.txt.gz.partac
)
if [[ ! -f "$GNOMAD_OUT" ]] || [[ $(stat -c%s "$GNOMAD_OUT") -lt 1000 ]]; then
    echo "  → Assembling gnomAD..."
    cat "${GNOMAD_PARTS[@]}" > "$GNOMAD_OUT"
    echo "  ✓ $GNOMAD_OUT"
else
    echo "  ✓ $GNOMAD_OUT (already assembled)"
fi

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

# ---------------------------------------------------------------------------
# 4. Verify resources
# ---------------------------------------------------------------------------
echo ""
echo "[4/4] Verifying resources..."

REQUIRED_FILES=(
    "resources/hg38/Homo_sapiens_assembly38.fasta"
    "resources/hg38/Homo_sapiens_assembly38.fasta.bwt"
    "resources/hg38/Homo_sapiens_assembly38.fasta.fai"
    "resources/hg38/Homo_sapiens_assembly38.dict"
    "resources/hg38/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
    "resources/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    "resources/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    "resources/hg38/1KG.20160622.strict_mask.hg38_GRCh38.fa.gz"
    "resources/hg38/gnomAD.r2.1.1.AFover0.001.snps.txt.gz"
    "resources/hg38/PON.q20q20.05.5.fa"
    "resources/hg38/PON.q20q20.05.5.fa.fai"
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
    echo "  Run sync_resources.sh from your LOCAL machine first:"
    echo "  bash scripts/sync_resources.sh <SERVER_USER>@<SERVER_HOST>:<SERVER_PATH>"
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
