#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-caddsv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

include { CADDSV } from './subworkflows/local/caddsv/main'

workflow {

    /*
     * Validate required params
     */
    if( !params.input )           error "Missing required parameter: --input"
    if( !params.annotations_dir ) error "Missing required parameter: --annotations_dir"
    if( !params.models_dir )      error "Missing required parameter: --models_dir"

    /*
     * Input channel:
     *  - If CSV: expects columns `sample` and `bed`
     *  - Else: assume it is a BED file path
     */
    ch_sv_bed = params.input.toString().endsWith('.csv')
        ? Channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def meta = [id: row.sample]
                def bed = file(row.bed, checkIfExists: true)
                tuple(meta, bed)
            }
        : Channel
            .fromPath(params.input, checkIfExists: true)
            .map { bed ->
                def meta = [id: bed.baseName.replaceAll(/\\.bed\$/, '')]
                tuple(meta, bed)
            }

    // Reference channels
    ch_annotations = Channel.fromPath(params.annotations_dir, checkIfExists: true, type: 'dir').first()
    ch_models      = Channel.fromPath(params.models_dir,      checkIfExists: true, type: 'dir').first()
    ch_scripts     = Channel.fromPath("${projectDir}/bin",    checkIfExists: true, type: 'dir').first()
    ch_genome      = Channel.fromPath("${projectDir}/assets/caddsv/hg38.genome", checkIfExists: true).first()

    // annotation header expected by score step
    ch_header = params.annotations_dir 
        ? Channel.fromPath("${params.annotations_dir}/header.txt", checkIfExists: true).first()
        : Channel.fromPath("${projectDir}/assets/caddsv/annotation_header.txt", checkIfExists: true).first()

    CADDSV(
        ch_sv_bed,
        ch_annotations,
        ch_models,
        ch_scripts,
        ch_genome
    )
}