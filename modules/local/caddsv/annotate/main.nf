/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CADDSV_ANNOTATE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CADDSV_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/53/5350a3386805dc54762f0bfd98ff3cc029ef50e9a536ad6f485f61d2564868c5/data' :
        'community.wave.seqera.io/library/bedops_bedtools_htslib_samtools_pruned:57bad587dad04e53' }"

    input:
    tuple val(meta),
          path(bed_wchr,      stageAs: 'bed_wchr.bed'),
          path(bed_nochr,     stageAs: 'bed_nochr.bed'),
          path(merged_wchr,   stageAs: 'merged_wchr.bed'),
          path(merged_nochr,  stageAs: 'merged_nochr.bed')
    path annotations_dir, stageAs: 'annotations'
    path scripts_dir,     stageAs: 'scripts'

    output:
    tuple val(meta), path("*_annotated.tsv"), emit: annotated
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env bash
    set -euo pipefail

    prefix="${prefix}"

    BED_WCHR="bed_wchr.bed"
    BED_NOCHR="bed_nochr.bed"
    MERG_WCHR="merged_wchr.bed"
    MERG_NOCHR="merged_nochr.bed"

    nlines=\$(wc -l < "\$BED_WCHR" | tr -d ' ')
    if [[ "\$nlines" -eq 0 ]]; then
      echo "ERROR: empty input BED: \$BED_WCHR" >&2
      exit 1
    fi

    # Helper functions
    dots_1col() { awk -v n="\$1" 'BEGIN{for(i=1;i<=n;i++)print "."}'; }
    dots_kcols() { awk -v n="\$1" -v k="\$2" 'BEGIN{OFS="\\t"; for(i=1;i<=n;i++){for(j=1;j<=k;j++){printf "%s%s",".",(j<k?OFS:ORS)}}}'; }

    # ========== CADD_PC_PhyloP (Snakemake: cadd_PC_phylop) ==========
    (while read -r line; do
        chr=\$(echo "\$line" | awk '{print \$1}')
        start=\$(echo "\$line" | awk '{print \$2}')
        end=\$(echo "\$line" | awk '{print \$3}')
        region="\${chr}:\${start}-\$((end+1))"
        bedtools map \\
          -b <(tabix "annotations/PhastCons/CADD_PC_PhyloP_scores.bed.gz" "\$region" 2>/dev/null | cat annotations/dummy8.bed - ) \\
          -a <(echo "\$line" | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$3}') \\
          -c 4,4,5,5,6,6,7,7,8,8 -o max,sum,max,sum,max,sum,max,sum,max,sum
      done < "\$BED_WCHR") > "\${prefix}_CADD_PC_PhyloP_maxsum.bed"

    # ========== CADD2 Count (Snakemake: cadd2) ==========
    (while read -r line; do
        chr=\$(echo "\$line" | awk '{print \$1}')
        start=\$(echo "\$line" | awk '{print \$2}')
        end=\$(echo "\$line" | awk '{print \$3}')
        region="\${chr}:\${start}-\${end}"
        paste <(echo "\$line" | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$3}') <(tabix "annotations/CADD/CADD_GRCh38-v1.5.bedGraph_90q_12.bed.gz" "\$region" 2>/dev/null | wc -l)
    done < "\$BED_WCHR") > "\${prefix}_cadd2_count.bed"

    # ========== CCR (Snakemake: ccr) ==========
    bedtools map -a "\$BED_NOCHR" -b "annotations/ccr/ccrs.all.bed.gz" -c 4 -o max > "\${prefix}_ccr_mean.bed"

    # ========== ChromHMM (Snakemake: chromHMM_MAX) ==========
    (while read -r line; do
        chr=\$(echo "\$line" | awk '{print \$1}')
        start=\$(echo "\$line" | awk '{print \$2}')
        end=\$(echo "\$line" | awk '{print \$3}')
        region="\${chr}:\${start}-\$((end+1))"
        bedtools map \\
          -b <(cat annotations/dummy_chrhmm.bed; tabix "annotations/chromhmm/chromHMM_GRCh38.bg.gz" "\$region" 2>/dev/null ) \\
          -a <(echo "\$line" | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$3}') \\
          -c 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28 -o max
      done < "\$BED_NOCHR") > "\${prefix}_chromHMM_max.bed"

    # ========== CTCF (Snakemake: CTCF) ==========
    bedtools coverage -b "annotations/CTCF/Wangetal_hg38.bed" -a "\$BED_WCHR" > "\${prefix}_ctcf.bed"

    # ========== Genomegitar DI (Snakemake: genomegitar1, genomegitar2, genomegitar3) ==========
    GG_DATASETS="GSM1055800_DI GSM1055805_DI GSM1081530_DI GSM1267196_DI GSM1267200_DI GSM1294038_DI GSM1294039_DI GSM1551599_DI GSM1551629_DI GSM1608505_DI GSM1718021_DI GSM1906332_DI GSM1906333_DI GSM1906334_DI GSM1909121_DI GSM455133_DI GSM862723_DI GSM862724_DI GSM927075_DI"
    
    for gg in \$GG_DATASETS; do
      bedtools map -a "\$BED_WCHR" -b "annotations/genomegitar/\${gg}/DI_sort.bed" -c 4,4 -o max,min > "\${prefix}_genomegitar_\${gg}.bed" 2>/dev/null || \\
        awk '{print \$0"\\t.\\t."}' "\$BED_WCHR" > "\${prefix}_genomegitar_\${gg}.bed"
    done

    # Paste all genomegitar files (Snakemake: genomegitar2)
    paste \${prefix}_genomegitar_*.bed > "\${prefix}_DI.bed"

    # Extract min/max (Snakemake: genomegitar3) - EXACT awk logic
    cut -f4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60,64,65,69,70,74,75,79,80,84,85,89,90 "\${prefix}_DI.bed" \\
      | awk '{m=\$1;for(i=1;i<=NF;i++)if(\$i<m)m=\$i;print m}' > "\${prefix}_DI_min.bed"
    cut -f4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60,64,65,69,70,74,75,79,80,84,85,89,90 "\${prefix}_DI.bed" \\
      | awk '{m=\$1;for(i=1;i<=NF;i++)if(\$i>m)m=\$i;print m}' > "\${prefix}_DI_max.bed"

    # ========== ENCODE (Snakemake: encode, encode2) ==========
    ENCODE_MARKS="DNase-seq H2AFZ H3K27ac H3K27me3 H3K36me3 H3K4me1 H3K4me2 H3K4me3 H3K79me2 H3K9ac H3K9me3 H4K20me1 totalRNA-seq"
    for mark in \$ENCODE_MARKS; do
      tabix "annotations/encode/\${mark}/\${mark}_merged_90quant.bed.gz" -R "\$MERG_NOCHR" 2>/dev/null \\
        | cat annotations/dummy5_nochr.bed - \\
        | bedtools map -a "\$BED_NOCHR" -b - -c 4,4 -o max,sum > "\${prefix}_encode_\${mark}_mean.bed"
    done

    # Paste all encode files (Snakemake: encode2)
    paste "\${prefix}_encode_DNase-seq_mean.bed" \\
          "\${prefix}_encode_H2AFZ_mean.bed" \\
          "\${prefix}_encode_H3K27ac_mean.bed" \\
          "\${prefix}_encode_H3K27me3_mean.bed" \\
          "\${prefix}_encode_H3K36me3_mean.bed" \\
          "\${prefix}_encode_H3K4me1_mean.bed" \\
          "\${prefix}_encode_H3K4me2_mean.bed" \\
          "\${prefix}_encode_H3K4me3_mean.bed" \\
          "\${prefix}_encode_H3K79me2_mean.bed" \\
          "\${prefix}_encode_H3K9ac_mean.bed" \\
          "\${prefix}_encode_H3K9me3_mean.bed" \\
          "\${prefix}_encode_H4K20me1_mean.bed" \\
          "\${prefix}_encode_totalRNA-seq_mean.bed" > "\${prefix}_encode.bed"

    # ========== EP (Snakemake: EP) ==========
    bedtools closest -d -t first -a "\$BED_WCHR" -b "annotations/enhancer-promoter-links/sorted_encode.bed" \\
      | Rscript --vanilla "scripts/annotateHIC.R" stdin "\${prefix}_EP.bed"

    # ========== FIRE (Snakemake: fire, fire2) ==========
    for cellline in gm12878 msc mes imr90 h1; do
      bedtools map -a "\$BED_NOCHR" -b "annotations/FIRE/fire_\${cellline}.bed" -c 4,4 -o max,min > "\${prefix}_fire_\${cellline}.bed"
    done

    # Paste fire files (Snakemake: fire2)
    paste "\${prefix}_fire_gm12878.bed" \\
          "\${prefix}_fire_msc.bed" \\
          "\${prefix}_fire_mes.bed" \\
          "\${prefix}_fire_imr90.bed" \\
          "\${prefix}_fire_h1.bed" > "\${prefix}_fire.bed"

    # ========== GC (Snakemake: GC) ==========
    (while read -r line; do
        chr=\$(echo "\$line" | awk '{print \$1}')
        start=\$(echo "\$line" | awk '{print \$2}')
        end=\$(echo "\$line" | awk '{print \$3}')
        region="\${chr}:\${start}-\$((end+1))"
        bedtools map \\
          -b <(tabix "annotations/GC/gc5Base.bedGraph.gz" "\$region" 2>/dev/null | cat annotations/dummy4.bed - ) \\
          -a <(echo "\$line" | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$3}') \\
          -c 4 -o mean
      done < "\$BED_WCHR") > "\${prefix}_gc.bed"

    # ========== Gene Model Tmp (Snakemake: gene_model_tmp) ==========
    tabix "annotations/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.gtf.gz" -R "\$MERG_NOCHR" \\
      | awk '{ if (\$0 ~ "transcript_id") print \$0; else print \$0" transcript_id \\"\\";"; }' \\
      | gtf2bed - \\
      | cut -f1,2,3,8,10 > "\${prefix}_gm_tmp.bed"

    # ========== Gene Model (Snakemake: gene_model) - NO -w flag! ==========
    paste \\
      <(bedtools coverage -b <(grep "exon" "\${prefix}_gm_tmp.bed") -a "\$BED_NOCHR" | cut -f1,2,3,5) \\
      <(bedtools coverage -b <(grep "transcript" "\${prefix}_gm_tmp.bed") -a "\$BED_NOCHR" | cut -f5) \\
      <(bedtools coverage -b <(grep "gene" "\${prefix}_gm_tmp.bed") -a "\$BED_NOCHR" | cut -f5) \\
      <(bedtools coverage -b <(grep "start_codon" "\${prefix}_gm_tmp.bed") -a "\$BED_NOCHR" | cut -f5) \\
      <(bedtools coverage -b <(grep "stop_codon" "\${prefix}_gm_tmp.bed") -a "\$BED_NOCHR" | cut -f5) \\
      <(bedtools coverage -b <(grep "three_prime_utr" "\${prefix}_gm_tmp.bed") -a "\$BED_NOCHR" | cut -f5) \\
      <(bedtools coverage -b <(grep "five_prime_utr" "\${prefix}_gm_tmp.bed") -a "\$BED_NOCHR" | cut -f5) \\
      <(bedtools coverage -b <(grep "CDS" "\${prefix}_gm_tmp.bed") -a "\$BED_NOCHR" | cut -f5) \\
      > "\${prefix}_genemodel.bed"

    # ========== Gene Model Dist (Snakemake: gene_model_dist) - NO -w flag! ==========
    grep 'exon' "annotations/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.bed" \\
      | bedtools closest -d -t first -b stdin -a "\$BED_NOCHR" | cut -f1,2,3,9 \\
      | paste - <(grep 'gene' "annotations/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.bed" | bedtools closest -d -t first -b stdin -a "\$BED_NOCHR" | cut -f9) \\
      | paste - <(grep 'start_codon' "annotations/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.bed" | bedtools closest -d -t first -b stdin -a "\$BED_NOCHR" | cut -f9) \\
      > "\${prefix}_genetic_dist.bed"

    # ========== Gene Names (Snakemake: gene_names) ==========
    bedtools map -b "annotations/ensembl_gff3/Homo_sapiens.GRCh38.96.chr.GENES.bed" -a "\$BED_NOCHR" -c 4 -o distinct > "\${prefix}_genenames.bed"

    # ========== PLI (Snakemake: pli) ==========
    Rscript --vanilla "scripts/PLIextract.R" "\${prefix}_genenames.bed" "annotations/gnomad/pli_exac.csv" "\${prefix}_pli.bed"

    # ========== GERP (Snakemake: gerp) ==========
    (while read -r line; do
        chr=\$(echo "\$line" | awk '{print \$1}')
        start=\$(echo "\$line" | awk '{print \$2}')
        end=\$(echo "\$line" | awk '{print \$3}')
        region="\${chr}:\${start}-\$((end+1))"
        bedtools map \\
          -b <(tabix "annotations/gerp/gerp_score2_hg38_MAM_90q.bed.gz" "\$region" 2>/dev/null | cat annotations/dummy5_nochr.bed - ) \\
          -a <(echo "\$line" | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$3}') \\
          -c 4 -o max
      done < "\$BED_NOCHR") > "\${prefix}_gerp_mean.bed"

    # ========== GERP2 Count (Snakemake: gerp2) ==========
    (while read -r line; do
        chr=\$(echo "\$line" | awk '{print \$1}')
        start=\$(echo "\$line" | awk '{print \$2}')
        end=\$(echo "\$line" | awk '{print \$3}')
        region="\${chr}:\${start}-\${end}"
        paste <(echo "\$line" | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$3}') <(tabix "annotations/gerp/gerp_score2_hg38_MAM_90q.bed.gz" "\$region" 2>/dev/null | wc -l)
    done < "\$BED_NOCHR") > "\${prefix}_gerp2_count.bed"

    # ========== HIC Encode (Snakemake: HIC_encode, mergetad) - stdin not /dev/stdin ==========
    for cell in A549 Caki2; do
      for ttype in nested tad; do
        bedtools closest -d -t first -a "\$BED_WCHR" -b "annotations/Encode-HIC/\${cell}/sorted_\${ttype}.bed.gz" \\
          | Rscript --vanilla "scripts/annotateHIC.R" stdin "\${prefix}_encode_\${cell}_\${ttype}_hic.bed"
      done
    done

    # Paste HIC encode files (Snakemake: mergetad)
    paste "\${prefix}_encode_A549_nested_hic.bed" \\
          "\${prefix}_encode_A549_tad_hic.bed" \\
          "\${prefix}_encode_Caki2_nested_hic.bed" \\
          "\${prefix}_encode_Caki2_tad_hic.bed" > "\${prefix}_HIC.bed"

    # ========== HIC_hESC (Snakemake: HIC_hESC) - stdin not /dev/stdin ==========
    bedtools closest -d -t first -a "\$BED_WCHR" -b "annotations/hic/hESC/combined/sorted.total.combined.domain" \\
      | Rscript --vanilla "scripts/annotateHIC.R" stdin "\${prefix}_HIC_hESC.bed"

    # ========== Microsynteny (Snakemake: microsynteny) - stdin not /dev/stdin ==========
    bedtools closest -d -t first -a "\$BED_NOCHR" -b "annotations/synteny/microsynteny.bed" \\
      | Rscript --vanilla "scripts/annotateHIC.R" stdin "\${prefix}_microsynteny.bed"

    # ========== MPC (Snakemake: MPC) ==========
    bedtools map -a "\$BED_NOCHR" -b "annotations/MPC/transcript_constraints_hg38liftover.bg.gz" -c 4 -o mean > "\${prefix}_MPC_mean.bed"

    # ========== RemapTF (Snakemake: RemapTF) ==========
    (while read -r line; do
        chr=\$(echo "\$line" | awk '{print \$1}')
        start=\$(echo "\$line" | awk '{print \$2}')
        end=\$(echo "\$line" | awk '{print \$3}')
        region="\${chr}:\${start}-\$((end+1))"
        bedtools map \\
          -b <(tabix "annotations/ReMap/ReMap2_overlapTF_hg38.bg.gz" "\$region" 2>/dev/null | cat annotations/dummy4_nochr.bed - ) \\
          -a <(echo "\$line" | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$3}') \\
          -c 4 -o mean
      done < "\$BED_NOCHR") > "\${prefix}_remapTF_mean.bed"

    # ========== Fantom5 (Snakemake: Fantom5_counts) ==========
    bedtools coverage -a "\$BED_WCHR" -b "annotations/fantom5/F5.hg38.enhancers_sort.bed" -counts > "\${prefix}_f5_counts.bed"

    # ========== DDD HI (Snakemake: HI) ==========
    bedtools map -a "\$BED_WCHR" -b "annotations/DDD_HI/hg38_HI_Predictions_version3_sort.bed" -c 5 -o max > "\${prefix}_dddhi.bed"

    # ========== DeepC (Snakemake: deepc) ==========
    (while read -r line; do
        chr=\$(echo "\$line" | awk '{print \$1}')
        start=\$(echo "\$line" | awk '{print \$2}')
        end=\$(echo "\$line" | awk '{print \$3}')
        region="\${chr}:\${start}-\${end}"
        paste <(echo "\$line" | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$3}') \\
          <(tabix "annotations/deepc/saliencies_merged_gm12878_5kb_10bp.bed.gz" "\$region" 2>/dev/null \\
            | awk 'BEGIN{maxVal="."}{if((maxVal==".")||(\$4>maxVal)){maxVal=\$4}}END{print maxVal}')
    done < "\$BED_WCHR") > "\${prefix}_deepc.bed"

    # ========== Ultraconserved (Snakemake: ultraconserved) ==========
    bedtools coverage -b "annotations/ultraconserved/ultraconserved_hg38_muliz120M_sort.bed" -a "\$BED_WCHR" > "\${prefix}_ultraconserved.bed"

    # ========== LINSIGHT (Snakemake: LINSIGHT) ==========
    (while read -r line; do
        chr=\$(echo "\$line" | awk '{print \$1}')
        start=\$(echo "\$line" | awk '{print \$2}')
        end=\$(echo "\$line" | awk '{print \$3}')
        region="\${chr}:\${start}-\${end}"
        paste <(echo "\$line" | awk 'BEGIN{OFS="\\t"}{print \$1,\$2,\$3}') \\
          <(tabix "annotations/linsight/LINSIGHT_hg38_sort.bed.gz" "\$region" 2>/dev/null \\
            | awk '{total+=\$4}END{printf("%f",total)}')
    done < "\$BED_WCHR") > "\${prefix}_linsight_sum.bed"

    # ========== FINAL ASSEMBLY (Snakemake: complete_script) - EXACT MATCH ==========
    paste \\
      <(cut -f1-11 "\${prefix}_CADD_PC_PhyloP_maxsum.bed") \\
      <(cut -f4 "\${prefix}_cadd2_count.bed") \\
      <(cut -f4 "\${prefix}_ccr_mean.bed") \\
      <(cut -f4-28 "\${prefix}_chromHMM_max.bed") \\
      <(cut -f4,5,7 "\${prefix}_ctcf.bed") \\
      <(cut -f1 "\${prefix}_DI_min.bed") \\
      <(cut -f1 "\${prefix}_DI_max.bed") \\
      <(cut -f4,5,9,10,14,15,19,20,24,25,29,30,34,35,39,40,44,45,49,50,54,55,59,60,64,65 "\${prefix}_encode.bed") \\
      <(cut -f4-7 "\${prefix}_EP.bed") \\
      <(cut -f4,5,9,10,14,15,19,20,24,25 "\${prefix}_fire.bed") \\
      <(cut -f4 "\${prefix}_gc.bed") \\
      <(cut -f4-11 "\${prefix}_genemodel.bed") \\
      <(cut -f4 "\${prefix}_gerp_mean.bed") \\
      <(cut -f4 "\${prefix}_gerp2_count.bed") \\
      <(cut -f4,5,6,7,11,12,13,14,18,19,20,21,25,26,27,28 "\${prefix}_HIC.bed") \\
      <(cut -f4-7 "\${prefix}_HIC_hESC.bed") \\
      <(cut -f4-7 "\${prefix}_microsynteny.bed") \\
      <(cut -f4 "\${prefix}_MPC_mean.bed") \\
      <(cut -f4 "\${prefix}_pli.bed") \\
      <(cut -f4,5,6 "\${prefix}_genetic_dist.bed") \\
      <(cut -f4 "\${prefix}_remapTF_mean.bed") \\
      <(cut -f4 "\${prefix}_f5_counts.bed") \\
      <(cut -f4 "\${prefix}_dddhi.bed") \\
      <(cut -f4 "\${prefix}_deepc.bed") \\
      <(cut -f4,5,7 "\${prefix}_ultraconserved.bed") \\
      <(cut -f4 "\${prefix}_linsight_sum.bed") \\
      | cat annotations/header.txt - > "\${prefix}_annotated.tsv"

    cat <<-END_VERSIONS > versions.yml
    "CADDSV_ANNOTATE":
      bedtools: "\$(bedtools --version | sed 's/bedtools v//')"
    END_VERSIONS
    """
}