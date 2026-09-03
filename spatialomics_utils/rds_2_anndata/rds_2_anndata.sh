#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 2 ]]; then
	echo "Usage: $(basename "$0") YOUR_RDS_PATH OUT_H5AD_PATH"
	exit 1
fi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
INPUT_RDS=$1
OUTPUT_H5AD=$2
OUTPUT_DIR=$(mktemp -d)

echo "Preparing temporary workspace: $OUTPUT_DIR"
trap 'rm -rf "$OUTPUT_DIR"' EXIT

echo "Splitting RDS into intermediate files..."
Rscript "$SCRIPT_DIR/break_rds.R" "$INPUT_RDS" "$OUTPUT_DIR"

echo "Building AnnData output: $OUTPUT_H5AD"
python3 "$SCRIPT_DIR/build_anndata.py" "$OUTPUT_DIR" "$OUTPUT_H5AD"

echo "Done."