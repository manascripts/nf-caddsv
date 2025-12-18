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
        }"

    input:
    tuple val(meta), path(bed)
    path genome

    output:
    tuple val(meta), path("${prefix}_wchr.bed")           , emit: bed_wchr
    tuple val(meta), path("${prefix}_nochr.bed")          , emit: bed_nochr
    tuple val(meta), path("${prefix}_merged_wchr.bed")    , emit: merged_wchr
    tuple val(meta), path("${prefix}_merged_nochr.bed")   , emit: merged_nochr
    tuple val(meta), path("${prefix}_100bpup_wchr.bed")   , emit: flank_up_wchr
    tuple val(meta), path("${prefix}_100bpup_nochr.bed")  , emit: flank_up_nochr
    tuple val(meta), path("${prefix}_100bpdown_wchr.bed") , emit: flank_down_wchr
    tuple val(meta), path("${prefix}_100bpdown_nochr.bed"), emit: flank_down_nochr
    path "versions.yml"                                   , emit: versions

    script:
    prefix = meta.id
    """
    #!/bin/bash
    set -euo pipefail

    # Prepare input: ensure we have chrom, start, end, svtype
    
    # Sort and normalize input
    awk 'BEGIN{OFS="\\t"} {
        chr = \$1
        start = \$2
        end = \$3
        svtype = (\$4 != "" ? \$4 : "DEL")
        
        # Ensure start < end
        if (start > end) { tmp = start; start = end; end = tmp }
        
        print chr, start, end, svtype
    }' ${bed} | sort -k1,1 -k2,2n > input_sorted.bed

    # Create with-chr version (add chr prefix if missing)
    awk 'BEGIN{OFS="\\t"} {
        chr = \$1
        if (chr !~ /^chr/) chr = "chr" chr
        print chr, \$2, \$3, \$4
    }' input_sorted.bed > ${prefix}_wchr.bed

    # Create no-chr version (remove chr prefix if present)
    awk 'BEGIN{OFS="\\t"} {
        chr = \$1
        sub(/^chr/, "", chr)
        print chr, \$2, \$3, \$4
    }' input_sorted.bed > ${prefix}_nochr.bed

    # Create merged versions (for annotations that need non-overlapping intervals)
    cut -f1-3 ${prefix}_wchr.bed | sort -k1,1 -k2,2n | bedtools merge -i - > ${prefix}_merged_wchr.bed
    cut -f1-3 ${prefix}_nochr.bed | sort -k1,1 -k2,2n | bedtools merge -i - > ${prefix}_merged_nochr.bed

    awk 'BEGIN{OFS="\\t"}{
        chr=\$1; len=\$2;
        if (chr ~ /^chr/) { print chr, len }
        else              { print "chr"chr, len }
    }' ${genome} | sort -k1,1 -k2,2n | uniq > genome_wchr_sorted.txt

    # Create 100bp upstream flanking regions
    bedtools flank -i ${prefix}_wchr.bed -g genome_wchr_sorted.txt -l 100 -r 0 2>/dev/null | \
        awk 'BEGIN{OFS="\\t"} \$3 > \$2 {print}' > ${prefix}_100bpup_wchr.bed || : > ${prefix}_100bpup_wchr.bed

    bedtools flank -i ${prefix}_nochr.bed -g ${genome} -l 100 -r 0 2>/dev/null | \
        awk 'BEGIN{OFS="\\t"} \$3 > \$2 {print}' > ${prefix}_100bpup_nochr.bed || : > ${prefix}_100bpup_nochr.bed

    # Create 100bp downstream flanking regions
    bedtools flank -i ${prefix}_wchr.bed -g genome_wchr_sorted.txt -l 0 -r 100 2>/dev/null | \
        awk 'BEGIN{OFS="\\t"} \$3 > \$2 {print}' > ${prefix}_100bpdown_wchr.bed || : > ${prefix}_100bpdown_wchr.bed

    bedtools flank -i ${prefix}_nochr.bed -g ${genome} -l 0 -r 100 2>/dev/null | \
        awk 'BEGIN{OFS="\\t"} \$3 > \$2 {print}' > ${prefix}_100bpdown_nochr.bed || : > ${prefix}_100bpdown_nochr.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
    END_VERSIONS
    """
}