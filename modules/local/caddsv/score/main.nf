
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CADDSV_SCORE - Score SVs using Random Forest models (Python)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CADDSV_SCORE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/scikit-learn:1.3.2' :
        'biocontainers/scikit-learn:1.3.2' }"

    input:
    tuple val(meta), path(annotated_main), path(annotated_up), path(annotated_down)
    path models_dir
    path scripts_dir
    path header_file

    output:
    tuple val(meta), path("*_scored.tsv")      , emit: scored
    tuple val(meta), path("*_scored_phred.tsv"), emit: scored_phred
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${scripts_dir}/scoring.py \\
        --annotated-main ${annotated_main} \\
        --annotated-up ${annotated_up} \\
        --annotated-down ${annotated_down} \\
        --models-dir ${models_dir} \\
        --header-file ${header_file} \\
        --output ${prefix}_scored.tsv \\
        --output-phred ${prefix}_scored_phred.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        scikit-learn: \$(python3 -c "import sklearn; print(sklearn.__version__)")
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
        joblib: \$(python3 -c "import joblib; print(joblib.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_scored.tsv
    touch ${prefix}_scored_phred.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.11.0
        scikit-learn: 1.3.2
        pandas: 2.0.0
        joblib: 1.3.0
    END_VERSIONS
    """
}
