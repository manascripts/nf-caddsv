/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CADDSV_ANNOTATE - Annotate SVs with genomic features
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs all CADD-SV annotation steps using bedtools and Python scripts.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CADDSV_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61571f6f07b3d-0' :
        'biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61571f6f07b3d-0' }"

    input:
    tuple val(meta), path(bed_wchr), path(bed_nochr), path(merged_wchr), path(merged_nochr)
    path annotations_dir
    path scripts_dir

    output:
    tuple val(meta), path("*_annotated.tsv"), emit: annotated
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash
    set -euo pipefail

    # Initialize output with coordinates
    cut -f1-4 ${bed_wchr} > ${prefix}_base.bed

    #==========================================================================
    # CONSERVATION & CONSTRAINT SCORES
    #==========================================================================

    echo "Annotating CTCF..."
    bedtools coverage -a ${bed_wchr} -b ${annotations_dir}/CTCF/Wangetal_hg38.bed 2>/dev/null | \\
        cut -f5 > ${prefix}_ctcf.bed || echo "." > ${prefix}_ctcf.bed

    echo "Annotating ultraconserved..."
    bedtools coverage -a ${bed_wchr} -b ${annotations_dir}/ultraconserved/ultraconserved_hg38_muliz120M_sort.bed 2>/dev/null | \\
        cut -f5 > ${prefix}_ultraconserved.bed || echo "." > ${prefix}_ultraconserved.bed

    echo "Annotating GC content..."
    tabix ${annotations_dir}/GC/gc5Base.bedGraph.gz -R ${merged_nochr} 2>/dev/null | \\
        bedtools map -a ${bed_nochr} -b stdin -c 4 -o mean 2>/dev/null | \\
        cut -f5 > ${prefix}_gc.bed || echo "." > ${prefix}_gc.bed

    echo "Annotating CADD/PhastCons/PhyloP..."
    tabix ${annotations_dir}/PhastCons/CADD_PC_PhyloP_scores.bed.gz -R ${merged_nochr} 2>/dev/null | \\
        bedtools map -a ${bed_nochr} -b stdin -c 4,5,6,7,8,9 -o max,max,max,sum,sum,sum 2>/dev/null | \\
        cut -f5-10 > ${prefix}_cadd_pc_phylop.bed || printf ".\\t.\\t.\\t.\\t.\\t.\\n" > ${prefix}_cadd_pc_phylop.bed

    echo "Annotating CADD counts..."
    tabix ${annotations_dir}/CADD/CADD_GRCh38-v1.5.bedGraph_90q_12.bed.gz -R ${merged_nochr} 2>/dev/null | \\
        bedtools coverage -a ${bed_nochr} -b stdin 2>/dev/null | \\
        cut -f5 > ${prefix}_cadd2_count.bed || echo "." > ${prefix}_cadd2_count.bed

    echo "Annotating GERP..."
    tabix ${annotations_dir}/gerp/gerp_score2_hg38_MAM_90q.bed.gz -R ${merged_nochr} 2>/dev/null | \\
        bedtools map -a ${bed_nochr} -b stdin -c 4 -o max 2>/dev/null | \\
        cut -f5 > ${prefix}_gerp_max.bed || echo "." > ${prefix}_gerp_max.bed

    tabix ${annotations_dir}/gerp/gerp_score2_hg38_MAM_90q.bed.gz -R ${merged_nochr} 2>/dev/null | \\
        bedtools coverage -a ${bed_nochr} -b stdin 2>/dev/null | \\
        cut -f5 > ${prefix}_gerp2_count.bed || echo "." > ${prefix}_gerp2_count.bed

    echo "Annotating LINSIGHT..."
    tabix ${annotations_dir}/linsight/LINSIGHT_hg38_sort.bed.gz -R ${merged_nochr} 2>/dev/null | \\
        bedtools map -a ${bed_nochr} -b stdin -c 4 -o sum 2>/dev/null | \\
        cut -f5 > ${prefix}_linsight.bed || echo "." > ${prefix}_linsight.bed

    echo "Annotating CCR..."
    tabix ${annotations_dir}/ccr/ccrs.all.bed.gz -R ${merged_nochr} 2>/dev/null | \\
        bedtools map -a ${bed_nochr} -b stdin -c 4 -o max 2>/dev/null | \\
        cut -f5 > ${prefix}_ccr.bed || echo "." > ${prefix}_ccr.bed

    echo "Annotating MPC..."
    tabix ${annotations_dir}/MPC/transcript_constraints_hg38liftover.bg.gz -R ${merged_nochr} 2>/dev/null | \\
        bedtools map -a ${bed_nochr} -b stdin -c 4 -o mean 2>/dev/null | \\
        cut -f5 > ${prefix}_mpc.bed || echo "." > ${prefix}_mpc.bed

    #==========================================================================
    # GENE MODEL ANNOTATIONS
    #==========================================================================

    echo "Annotating gene models..."
    tabix ${annotations_dir}/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.gtf.gz -R ${merged_wchr} 2>/dev/null > ${prefix}_gm_tmp.gtf || touch ${prefix}_gm_tmp.gtf

    for feature in exon transcript gene start_codon stop_codon three_prime_utr five_prime_utr CDS; do
        grep "\${feature}" ${prefix}_gm_tmp.gtf 2>/dev/null | \\
            bedtools coverage -a ${bed_wchr} -b stdin 2>/dev/null | \\
            cut -f5 > ${prefix}_gm_\${feature}.bed || echo "." > ${prefix}_gm_\${feature}.bed
    done

    echo "Annotating gene distance..."
    bedtools closest -a ${bed_wchr} -b ${annotations_dir}/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.bed -d 2>/dev/null | \\
        awk '{print \$NF}' > ${prefix}_gene_dist.bed || echo "." > ${prefix}_gene_dist.bed

    echo "Extracting gene names..."
    bedtools map -a ${bed_wchr} -b ${annotations_dir}/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.GENES.bed \\
        -c 4 -o distinct 2>/dev/null > ${prefix}_gene_names.bed || touch ${prefix}_gene_names.bed

    echo "Annotating pLI (using Python)..."
    python3 ${scripts_dir}/pli_extract.py \\
        --gene-names ${prefix}_gene_names.bed \\
        --pli-file ${annotations_dir}/gnomad/pli_exac.csv \\
        --output ${prefix}_pli.bed || echo "." > ${prefix}_pli.bed

    #==========================================================================
    # REGULATORY FEATURES
    #==========================================================================

    echo "Annotating ReMap TF..."
    tabix ${annotations_dir}/ReMap/ReMap2_overlapTF_hg38.bg.gz -R ${merged_nochr} 2>/dev/null | \\
        bedtools map -a ${bed_nochr} -b stdin -c 4 -o mean 2>/dev/null | \\
        cut -f5 > ${prefix}_remap.bed || echo "." > ${prefix}_remap.bed

    echo "Annotating FANTOM5..."
    bedtools coverage -a ${bed_wchr} -b ${annotations_dir}/fantom5/F5.hg38.enhancers_sort.bed 2>/dev/null | \\
        cut -f5 > ${prefix}_fantom5.bed || echo "." > ${prefix}_fantom5.bed

    echo "Annotating HI..."
    bedtools map -a ${bed_wchr} -b ${annotations_dir}/DDD_HI/hg38_HI_Predictions_version3_sort.bed \\
        -c 5 -o max 2>/dev/null | cut -f5 > ${prefix}_hi.bed || echo "." > ${prefix}_hi.bed

    echo "Annotating DeepC..."
    tabix ${annotations_dir}/deepc/saliencies_merged_gm12878_5kb_10bp.bed.gz -R ${merged_nochr} 2>/dev/null | \\
        bedtools map -a ${bed_nochr} -b stdin -c 4 -o max 2>/dev/null | \\
        cut -f5 > ${prefix}_deepc.bed || echo "." > ${prefix}_deepc.bed

    echo "Annotating chromHMM..."
    tabix ${annotations_dir}/chromhmm/chromHMM_GRCh38.bg.gz -R ${merged_wchr} 2>/dev/null > ${prefix}_chromhmm_tmp.bed || touch ${prefix}_chromhmm_tmp.bed
    if [ -s ${prefix}_chromhmm_tmp.bed ]; then
        for i in \$(seq 4 28); do
            bedtools map -a ${bed_wchr} -b ${prefix}_chromhmm_tmp.bed -c \$i -o max 2>/dev/null | \\
                cut -f5 > ${prefix}_chromhmm_\${i}.bed || echo "." > ${prefix}_chromhmm_\${i}.bed
        done
        paste ${prefix}_chromhmm_*.bed > ${prefix}_chromhmm.bed
    else
        # Create 25 empty columns
        yes "." | head -n \$(wc -l < ${bed_wchr}) | paste -d'\\t' \$(for i in \$(seq 1 25); do echo "-"; done) > ${prefix}_chromhmm.bed
    fi

    echo "Annotating ENCODE histone marks..."
    ENCODE_MARKS="DNase-seq H2AFZ H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me2 H3K4me3 H3K79me2 H3K9ac H3K9me3 H4K20me1 totalRNA-seq"
    for mark in \${ENCODE_MARKS}; do
        tabix ${annotations_dir}/encode/\${mark}/\${mark}_merged_90quant.bed.gz -R ${merged_nochr} 2>/dev/null | \\
            bedtools map -a ${bed_nochr} -b stdin -c 4,4 -o max,sum 2>/dev/null | \\
            cut -f5,6 > ${prefix}_encode_\${mark}.bed || printf ".\\t.\\n" > ${prefix}_encode_\${mark}.bed
    done
    paste ${prefix}_encode_*.bed > ${prefix}_encode.bed

    #==========================================================================
    # 3D GENOME FEATURES
    #==========================================================================

    echo "Annotating FIRE..."
    for cl in gm12878 msc mes imr90 h1; do
        bedtools map -a ${bed_wchr} -b ${annotations_dir}/FIRE/fire_\${cl}.bed -c 4,4 -o max,min 2>/dev/null | \\
            cut -f5,6 > ${prefix}_fire_\${cl}.bed || printf ".\\t.\\n" > ${prefix}_fire_\${cl}.bed
    done
    paste ${prefix}_fire_*.bed > ${prefix}_fire.bed

    echo "Annotating microsynteny..."
    bedtools closest -a ${bed_wchr} -b ${annotations_dir}/synteny/microsynteny.bed -d 2>/dev/null | \\
        awk '{print \$NF}' > ${prefix}_synteny.bed || echo "." > ${prefix}_synteny.bed

    echo "Annotating EP links..."
    bedtools closest -a ${bed_wchr} -b ${annotations_dir}/enhancer-promoter-links/sorted_encode.bed -d 2>/dev/null | \\
        awk '{print \$NF}' > ${prefix}_ep.bed || echo "." > ${prefix}_ep.bed

    #==========================================================================
    # COMBINE ALL ANNOTATIONS
    #==========================================================================

    echo "Combining all annotations..."
    paste \\
        ${prefix}_base.bed \\
        ${prefix}_ctcf.bed \\
        ${prefix}_ultraconserved.bed \\
        ${prefix}_gc.bed \\
        ${prefix}_cadd_pc_phylop.bed \\
        ${prefix}_cadd2_count.bed \\
        ${prefix}_gerp_max.bed \\
        ${prefix}_gerp2_count.bed \\
        ${prefix}_linsight.bed \\
        ${prefix}_ccr.bed \\
        ${prefix}_mpc.bed \\
        ${prefix}_gm_exon.bed \\
        ${prefix}_gm_transcript.bed \\
        ${prefix}_gm_gene.bed \\
        ${prefix}_gm_start_codon.bed \\
        ${prefix}_gm_stop_codon.bed \\
        ${prefix}_gm_three_prime_utr.bed \\
        ${prefix}_gm_five_prime_utr.bed \\
        ${prefix}_gm_CDS.bed \\
        ${prefix}_gene_dist.bed \\
        ${prefix}_pli.bed \\
        ${prefix}_remap.bed \\
        ${prefix}_fantom5.bed \\
        ${prefix}_hi.bed \\
        ${prefix}_deepc.bed \\
        ${prefix}_chromhmm.bed \\
        ${prefix}_encode.bed \\
        ${prefix}_fire.bed \\
        ${prefix}_synteny.bed \\
        ${prefix}_ep.bed \\
        > ${prefix}_annotated.tsv

    # Clean up intermediate files
    rm -f ${prefix}_*.bed ${prefix}_gm_tmp.gtf ${prefix}_chromhmm_tmp.bed 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
        tabix: \$(tabix --version 2>&1 | head -n1 | sed 's/tabix (htslib) //')
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed 's/bedtools v//')
        tabix: \$(tabix --version 2>&1 | head -n1 | sed 's/tabix (htslib) //')
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}