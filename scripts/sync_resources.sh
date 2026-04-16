#!/usr/bin/env bash
# =============================================================================
# sync_resources.sh — Transfer resources/ to the server
#
# Run this script on your LOCAL machine (WSL2).
#
# Usage:
#   bash scripts/sync_resources.sh <user>@<host>:<remote_project_dir>
#
# Example:
#   bash scripts/sync_resources.sh user@server.example.com:/home/user/BSMN_pipeline_ij
#
# What is transferred (~12 GB):
#   resources/hg38/  — reference genome, known sites, gnomAD, masks
#   downloads/       — split files (gnomAD, PON) to be assembled on server
#
# Containers are NOT transferred (pulled on server via server_setup.sh).
# =============================================================================
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: bash scripts/sync_resources.sh <user>@<host>:<remote_project_dir>"
    echo ""
    echo "Example:"
    echo "  bash scripts/sync_resources.sh user@server.example.com:/home/user/BSMN_pipeline_ij"
    exit 1
fi

REMOTE="$1"
LOCAL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "[sync] Source : $LOCAL_DIR/resources/ + $LOCAL_DIR/downloads/"
echo "[sync] Target : $REMOTE/resources/ + $REMOTE/downloads/"
echo "[sync] Starting rsync (this may take a while for ~12 GB)..."
echo ""

rsync -avhP \
    --stats \
    "$LOCAL_DIR/resources/" \
    "$REMOTE/resources/"

echo ""
echo "[sync] Transferring split resource files (downloads/)..."
rsync -avhP \
    --stats \
    "$LOCAL_DIR/downloads/" \
    "$REMOTE/downloads/"

echo ""
echo "[sync] Done."
echo ""
echo "Next step — run on the server:"
echo "  cd <remote_project_dir>"
echo "  bash scripts/server_setup.sh"
