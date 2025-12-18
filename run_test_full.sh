#!/bin/bash
set -euo pipefail

export PATH=$PATH:/omics/odcf/analysis/OE0415_projects/nb_lrs/LRS/2025-02-01_NB_LRS/01-src/nf

cd "$(dirname "$0")"

# Verify references exist
if [ ! -d "caddsv_references/annotations" ]; then
    echo "ERROR: Annotations not found. Run ./run_test_download.sh first"
    exit 1
fi

if [ ! -d "caddsv_references/models" ]; then
    echo "ERROR: Models not found. Run ./run_test_download.sh first"
    exit 1
fi

nextflow run main.nf \
    -profile singularity \
    --input test/data/test_svs.bed \
    --annotations_dir ./caddsv_references/annotations \
    --models_dir ./caddsv_references/models \
    --outdir test/results \
    -resume

