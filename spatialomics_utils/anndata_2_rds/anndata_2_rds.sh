#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 2 ]]; then
	echo "Usage: $(basename "$0") YOUR_ADATA_PATH OUT_RDS_PATH"
	exit 1
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ADATA_PATH=$1
OUTPUT_RDS=$2
OUTPUT_DIR=$(mktemp -d)

echo "Preparing temporary workspace: $OUTPUT_DIR"
trap 'rm -rf "$OUTPUT_DIR"' EXIT

echo "Splitting AnnData into intermediate files..."
python3 "$SCRIPT_DIR/break_anndata.py" "$ADATA_PATH" "$OUTPUT_DIR"

echo "Building RDS output: $OUTPUT_RDS"
Rscript "$SCRIPT_DIR/build_rds.R" "$OUTPUT_DIR" "$OUTPUT_RDS"

echo "Done."