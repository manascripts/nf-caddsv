#!/bin/bash
set -euo pipefail

export PATH=$PATH:/omics/odcf/analysis/OE0415_projects/nb_lrs/LRS/2025-02-01_NB_LRS/01-src/nf

cd "$(dirname "$0")"

nextflow run test_download.nf \
    -profile singularity \
    --caddsv_cache_dir ./results/caddsv_references \
    -resume
