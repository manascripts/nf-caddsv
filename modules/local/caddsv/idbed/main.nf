/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CADDSV_IDBED - Create input/id_<sample>.bed expected by the official CADD-SV scoring.R
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CADDSV_IDBED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7d/7df273d12f0c4d8539440b68876edf39b739cb78bb806418c5b5d057fe11bdbd/data' :
        'community.wave.seqera.io/library/bedtools:2.31.1--7c4ce4cb07c09ee4' }"
    // same conatiner as prep_bed
    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("id_${meta.id}.bed"), emit: idbed

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Ensure 5 columns: chr start end SVTYPE name
    awk 'BEGIN{OFS="\\t"}{
        chr=\$1; start=\$2; end=\$3;
        svt=(NF>=4 && \$4!="" ? \$4 : "DEL");
        name=(NF>=5 && \$5!="" ? \$5 : ".");
        if (start > end) { tmp=start; start=end; end=tmp }
        print chr, start, end, svt, name
    }' ${bed} > id_${meta.id}.bed
    """ 
}