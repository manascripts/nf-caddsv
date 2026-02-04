/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CADDSV_SCORE - Score SVs using Random Forest models (R, faithful to Snakemake)
    CRITICAL: Snakemake rule scoring (line 1067-1079) has counter-intuitive naming:
      input:
        span = matrix.bed
        flank_up = matrix_100bpdown.bed      # Uses DOWN matrix for UP context
        flank_down = matrix_100bpup.bed      # Uses UP matrix for DOWN context
    
    This is scientifically intentional - maintain exact Snakemake order.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CADDSV_SCORE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a3/a3518b248237cdea316a05907915a81fa26d6999c5cf03ad9fde4f8f1536e13b/data' :
        'community.wave.seqera.io/library/bedtools_r-base_r-data.table_r-optparse_pruned:c587fff83c849578' }"
    
    input:
    tuple val(meta), 
        path(original_bed), 
        path(matrix_main), 
        path(matrix_up), path(matrix_down)
    path models_dir
    path scripts_dir
    path genome_file
    path annotations_dir, stageAs: 'annotations'

    output:
    tuple val(meta), path("${meta.id}_score.bed"), emit: score
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def name = meta.id
    """
    mkdir -p input output "${name}"
    # scoring.R expects to find models in models/ subdirectory via relative path
    # It reads this internally, not via command line args
    # ln -s \$(realpath ${models_dir}) models

    # scoring.R reads: input/id_<name>.bed
    cp ${original_bed} "input/id_${name}.bed"

    # Matrix files must be in <name>/ directory with exact names
    cp ${matrix_main} "${name}/matrix.bed"
    cp ${matrix_up} "${name}/matrix_100bpup.bed"
    cp ${matrix_down} "${name}/matrix_100bpdown.bed"
    

    # Snakemake command (line 1078):
    # Rscript --vanilla scripts/scoring.R {params.name} {input.span} {input.flank_up} {input.flank_down} {input.genome} {output}
    # Args breakdown:
    #   args[1] = name (dataset identifier)
    #   args[2] = span matrix path
    #   args[3] = flank_up matrix path (actually downstream in Snakemake)
    #   args[4] = flank_down matrix path (actually upstream in Snakemake)
    #   args[5] = genome file path
    #   args[6] = output file path
    
    # scoring.R only uses args[1] (name) and args[5] (genome)
    # args[2-4] and args[6] are documented but the function reads from fixed paths
    Rscript --vanilla ${scripts_dir}/scoring.R \\
        "${name}" \\
        "${matrix_main}" \\
        "${matrix_down}" \\
        "${matrix_up}" \\
        "${genome_file}" \\
        "${name}.score"

    # ===== Snakemake rule sort =====
    # bedtools sort -i {input.score} | cat {input.header} - > {output}
    bedtools sort -i "output/${name}.score" \\
        | cat annotations/header_final.txt - \\
        > "${name}_score.bed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/R version //' | sed 's/ (.*//')
    END_VERSIONS
    """
}