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

echo "[sync] Source : $LOCAL_DIR/resources/"
echo "[sync] Target : $REMOTE/resources/"
echo "[sync] Starting rsync (this may take a while for ~12 GB)..."
echo ""

rsync -avhP \
    --stats \
    "$LOCAL_DIR/resources/" \
    "$REMOTE/resources/"

echo ""
echo "[sync] Done. Resources transferred to $REMOTE/resources/"
echo ""
echo "Next step — run on the server:"
echo "  cd <remote_project_dir>"
echo "  bash scripts/server_setup.sh"
