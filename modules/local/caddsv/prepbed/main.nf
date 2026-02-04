/*
 * CADDSV_PREPBED - Prepare BED files for CADD-SV annotation
 * 
 * Creates multiple BED file versions:
 * - With/without chr prefix (different annotations use different conventions)
 * - Merged intervals for overlapping SVs
 * - 100bp flanking regions (upstream and downstream)
 */

process CADDSV_PREPBED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7df273d12f0c4d8539440b68876edf39b739cb78bb806418c5b5d057fe11bdbd/data' :
        'community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4' }"

    input:
    tuple val(meta), path(bed)
    path genome

    output:
    tuple val(meta), path("${meta.id}_sorted_original.bed"), emit: sorted_original
    // Main SV region BED files
    tuple val(meta), 
          path("${meta.id}_wchr.bed"), 
          path("${meta.id}_nochr.bed"),
          path("${meta.id}_merged.bed"),
          path("${meta.id}_nochr_merged.bed"), emit: bed_main
    
    // 100bp upstream flank
    tuple val(meta),
          path("${meta.id}_100bpup.bed"),
          path("${meta.id}_100bpup_nchr.bed"),
          path("${meta.id}_100bpup_merged.bed"),
          path("${meta.id}_100bpup_nchr_merged.bed"), emit: bed_flank_up
    
    // 100bp downstream flank
    tuple val(meta),
          path("${meta.id}_100bpdown.bed"),
          path("${meta.id}_100bpdown_nchr.bed"),
          path("${meta.id}_100bpdown_merged.bed"),
          path("${meta.id}_100bpdown_nchr_merged.bed"), emit: bed_flank_down
    
    path "versions.yml", emit: versions

    script:
    prefix = meta.id
    """
    #!/bin/bash
    set -euo pipefail

    # ===== Sort the ORIGINAL input BED (all columns) and save for scoring =====
    bedtools sort -i ${bed} > "${prefix}_sorted_original.bed"

    # ===== prep_chr1 - Extract chr,start,end from SORTED original =====
    cut -f1,2,3 "${prefix}_sorted_original.bed" > "${prefix}_wchr.bed"

    # ===== prep_chr2 - Remove chr prefix =====
    sed 's/^chr//g' "${prefix}_wchr.bed" > "${prefix}_nochr.bed"

    # ===== prep_merg1 =====
    bedtools merge -i "${prefix}_nochr.bed" > "${prefix}_nochr_merged.bed"
    bedtools merge -i "${prefix}_wchr.bed" > "${prefix}_merged.bed"
    
    # ========== prep_chr_100bpup ==========
    awk 'BEGIN{OFS="\\t"}{if(\$2==0)\$2+=1; print \$0}' "${prefix}_wchr.bed" > "${prefix}_100bpup_tmp.bed"
    bedtools flank -i "${prefix}_100bpup_tmp.bed" -g ${genome} -l 100 -r 0 | bedtools sort > "${prefix}_100bpup.bed"
    sed 's/^chr//g' "${prefix}_100bpup.bed" > "${prefix}_100bpup_nchr.bed"

    # ========== prep_merg1_100bpup ==========
    bedtools merge -i "${prefix}_100bpup_nchr.bed" > "${prefix}_100bpup_nchr_merged.bed"
    bedtools merge -i "${prefix}_100bpup.bed" > "${prefix}_100bpup_merged.bed"

    # ========== prep_chr_100bpdown ==========
    bedtools flank -i "${prefix}_wchr.bed" -g ${genome} -l 0 -r 100 | bedtools sort > "${prefix}_100bpdown.bed"
    sed 's/^chr//g' "${prefix}_100bpdown.bed" > "${prefix}_100bpdown_nchr.bed"

    # ========== prep_merg1_100bpdown ==========
    bedtools merge -i "${prefix}_100bpdown_nchr.bed" > "${prefix}_100bpdown_nchr_merged.bed"
    bedtools merge -i "${prefix}_100bpdown.bed" > "${prefix}_100bpdown_merged.bed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}