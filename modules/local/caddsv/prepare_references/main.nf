/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CADDSV_PREPARE_REFERENCES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Downloads the CADD-SV dependencies tarball and models
    URL: https://kircherlab.bihealth.org/download/CADD-SV/v1.1/dependencies.tar.gz
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CADDSV_PREPARE_REFERENCES {
    tag "download_caddsv_refs"
    label 'process_low'
    storeDir "${params.caddsv_cache_dir ?: 'caddsv_references'}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7c/7c048e430f136e1910a1564516ef6837a6ed3c951c7e7bf8ff88a87ac8642cf1/data' :
        'community.wave.seqera.io/library/curl_tar:5d98681149abf1eb' }"

    output:
    path "annotations", emit: annotations_dir
    path "models"     , emit: models_dir
    path "versions.yml", emit: versions

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Detect tag (from env or default)
    TAG="\${CADDSV_TAG:-v1.1}"

    # Try requested tag
    curl -sSL "https://kircherlab.bihealth.org/download/CADD-SV/\$TAG/dependencies.tar.gz" -o dependencies.tar.gz || true

    # Fallback to v1.1 if download failed or file is empty
    if [[ ! -s dependencies.tar.gz ]]; then
        echo "Falling back to v1.1..."
        curl -sSL "https://kircherlab.bihealth.org/download/CADD-SV/v1.1/dependencies.tar.gz" -o dependencies.tar.gz || { echo "Download failed" >&2; exit 1; }
    fi

    # Extract annotations
    tar -xzf dependencies.tar.gz
    rm dependencies.tar.gz

    # Download models for the same tag
    mkdir -p models
    curl -s "https://api.github.com/repos/kircherlab/CADD-SV/contents/models?ref=\$TAG" \
      | grep '"download_url"' | cut -d '"' -f 4 \
      | while read url; do
          curl -sSL "\$url" -o "models/\$(basename "\$url")"
        done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tag: "\$TAG"
    END_VERSIONS
    """
}
