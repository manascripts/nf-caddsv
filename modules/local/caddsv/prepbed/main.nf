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
    tuple val(meta), path("${prefix}_wchr.bed"),         path("${prefix}_nochr.bed"),          path("${prefix}_merged_wchr.bed"),    path("${prefix}_merged_nochr.bed"),    emit: bed_main
    tuple val(meta), path("${prefix}_100bpup_wchr.bed"), path("${prefix}_100bpup_nochr.bed"),  path("${prefix}_100bpup_wchr.bed"),  path("${prefix}_100bpup_nochr.bed"),  emit: bed_flank_up
    tuple val(meta), path("${prefix}_100bpdown_wchr.bed"), path("${prefix}_100bpdown_nochr.bed"), path("${prefix}_100bpdown_wchr.bed"), path("${prefix}_100bpdown_nochr.bed"), emit: bed_flank_down
    path "versions.yml"                                   , emit: versions

    script:
    prefix = meta.id
    """
    #!/bin/bash
    set -euo pipefail

    # Sort and normalize input
    awk 'BEGIN{OFS="\\t"} {
        chr = \$1
        start = \$2
        end = \$3
        svtype = (\$4 != "" ? \$4 : "DEL")
        if (start > end) { tmp = start; start = end; end = tmp }
        print chr, start, end, svtype
    }' ${bed} | sort -k1,1 -k2,2n > input_sorted.bed

    # with-chr
    awk 'BEGIN{OFS="\\t"} {
        chr = \$1
        if (chr !~ /^chr/) chr = "chr" chr
        print chr, \$2, \$3, \$4
    }' input_sorted.bed > ${prefix}_wchr.bed

    # no-chr
    awk 'BEGIN{OFS="\\t"} {
        chr = \$1
        sub(/^chr/, "", chr)
        print chr, \$2, \$3, \$4
    }' input_sorted.bed > ${prefix}_nochr.bed

    # merged
    cut -f1-3 ${prefix}_wchr.bed | sort -k1,1 -k2,2n | bedtools merge -i - > ${prefix}_merged_wchr.bed
    cut -f1-3 ${prefix}_nochr.bed | sort -k1,1 -k2,2n | bedtools merge -i - > ${prefix}_merged_nochr.bed

    # -------------------------------------------------------------------
    # FIX: build chr-prefixed genome once, do NOT generate chrchr contigs
    # -------------------------------------------------------------------
    awk 'BEGIN{OFS="\\t"}{
        chr=\$1; len=\$2;
        if (chr ~ /^chr/) print chr, len;
        else              print "chr"chr, len;
    }' ${genome} | sort -k1,1 -k2,2n | uniq > genome_wchr_sorted.txt

    # flanks
    bedtools flank -i ${prefix}_wchr.bed   -g genome_wchr_sorted.txt -l 100 -r 0   2>/dev/null | awk 'BEGIN{OFS="\\t"} \$3>\$2{print}' > ${prefix}_100bpup_wchr.bed   || : > ${prefix}_100bpup_wchr.bed
    bedtools flank -i ${prefix}_nochr.bed  -g ${genome}              -l 100 -r 0   2>/dev/null | awk 'BEGIN{OFS="\\t"} \$3>\$2{print}' > ${prefix}_100bpup_nochr.bed  || : > ${prefix}_100bpup_nochr.bed
    bedtools flank -i ${prefix}_wchr.bed   -g genome_wchr_sorted.txt -l 0   -r 100 2>/dev/null | awk 'BEGIN{OFS="\\t"} \$3>\$2{print}' > ${prefix}_100bpdown_wchr.bed || : > ${prefix}_100bpdown_wchr.bed
    bedtools flank -i ${prefix}_nochr.bed  -g ${genome}              -l 0   -r 100 2>/dev/null | awk 'BEGIN{OFS="\\t"} \$3>\$2{print}' > ${prefix}_100bpdown_nochr.bed|| : > ${prefix}_100bpdown_nochr.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}