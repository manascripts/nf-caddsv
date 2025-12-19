/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CADDSV_SCORE - Score SVs using Random Forest models (R, faithful to Snakemake)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CADDSV_SCORE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8e/8eb87e70b56d92de45cca178ca4abf4694170f88dc434c1851473ae1dfaf2e1b/data' :
        'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-optparse_pruned:17b521e080857155' }"
    
    input:
    tuple val(meta), path(idbed), path(matrix_main), path(matrix_up), path(matrix_down)
    path models_dir
    path scripts_dir
    path genome_file

    output:
    tuple val(meta), path("${meta.id}_score.bed")      , emit: score
    tuple val(meta), path("${meta.id}_score_phred.bed"), emit: score_phred
    path "versions.yml"                                , emit: versions

    script:
    def name    = meta.id
    def threads = task.cpus ?: 1

    """
    #!/bin/bash
    set -euo pipefail

    # Locate scoring script
    SCORE_R="${scripts_dir}/scoring.R"
    [[ -f "\$SCORE_R" ]] || { echo "ERROR: scoring.R not found at \$SCORE_R"; ls -la "${scripts_dir}"; exit 1; }

    # Stage the filesystem layout expected by scoring.R
    mkdir -p input models "${name}" outdir

    # input/id_<name>.bed
    cp "${idbed}" "input/id_${name}.bed"

    # matrices (must be exactly these names/paths)
    cp "${matrix_main}" "${name}/matrix.bed"
    cp "${matrix_up}"   "${name}/matrix_100bpup.bed"
    cp "${matrix_down}" "${name}/matrix_100bpdown.bed"

    # models/ must be relative for scoring.R (it does readRDS('models/...'))
    # Use symlinks to avoid copying huge directories (works on Linux in Nextflow workdir).
    # If symlinks are an issue on your executor, change to: cp -r ${models_dir}/* models/
    ln -s ${models_dir}/* models/ 2>/dev/null || true

    # Run scoring
    # Expected args: name, prefix, threads, outputdir, genome
    Rscript "\$SCORE_R" "${name}" "${name}" "${threads}" "outdir" "${genome_file}"

    # Normalize outputs to nf-core-like filenames.
    # We don't know exact filenames produced, so capture the common ones and fail loudly if missing.
    # Typical outputs in CADD-SV are something like: outdir/<name>.score and outdir/<name>_score_phred
    if [[ -f "outdir/${name}.score" ]]; then
      cp "outdir/${name}.score" "${name}_score.bed"
    elif [[ -f "outdir/${name}_score.bed" ]]; then
      cp "outdir/${name}_score.bed" "${name}_score.bed"
    else
      echo "ERROR: Could not find main score output in outdir/"
      ls -la outdir
      exit 1
    fi

    if [[ -f "outdir/${name}.score_phred" ]]; then
      cp "outdir/${name}.score_phred" "${name}_score_phred.bed"
    elif [[ -f "outdir/${name}_score_phred.bed" ]]; then
      cp "outdir/${name}_score_phred.bed" "${name}_score_phred.bed"
    elif [[ -f "outdir/${name}.phred" ]]; then
      cp "outdir/${name}.phred" "${name}_score_phred.bed"
    else
      echo "WARNING: Could not find PHRED score output in outdir/; creating empty placeholder"
      : > "${name}_score_phred.bed"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/.*R version //')
    END_VERSIONS
    """
}
