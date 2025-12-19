/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CADDSV_ASSEMBLE
    Create CADD-SV matrix files expected by the official scoring.R:
      <sample>/matrix.bed
      <sample>/matrix_100bpup.bed
      <sample>/matrix_100bpdown.bed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CADDSV_ASSEMBLE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7df273d12f0c4d8539440b68876edf39b739cb78bb806418c5b5d057fe11bdbd/data' :
        'community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4' }"
    // same container as prepbed
    input:
    tuple val(meta), path(annotated_main), path(annotated_up), path(annotated_down)
    path header_file

    output:
    tuple val(meta), path("${meta.id}/matrix.bed")           , emit: matrix_main
    tuple val(meta), path("${meta.id}/matrix_100bpup.bed")   , emit: matrix_up
    tuple val(meta), path("${meta.id}/matrix_100bpdown.bed") , emit: matrix_down
    path "versions.yml"                                      , emit: versions

    script:
    """
    #!/bin/bash
    set -euo pipefail

    mkdir -p "${meta.id}"

    strip_header() {
      awk 'NR==1{
             # Drop a possible header line
             if (\$0 ~ /^#/ || tolower(\$1) ~ /chrom/ || tolower(\$1) == "chr") { next }
           }
           { print }
      ' "\$1"
    }

    # Official pipeline prepends header.txt to each matrix.
    cat "${header_file}" <(strip_header "${annotated_main}") > "${meta.id}/matrix.bed"
    cat "${header_file}" <(strip_header "${annotated_up}")   > "${meta.id}/matrix_100bpup.bed"
    cat "${header_file}" <(strip_header "${annotated_down}") > "${meta.id}/matrix_100bpdown.bed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: "builtin"
    END_VERSIONS
    """
}