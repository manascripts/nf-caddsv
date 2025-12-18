#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-caddsv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

include { CADDSV } from './subworkflows/local/caddsv/main'

workflow {

    log.info """
    ══════════════════════════════════════════════════════════════════════════
      nf-caddsv : CADD-SV Structural Variant Scoring Pipeline
    ══════════════════════════════════════════════════════════════════════════
    Input            : ${params.input}
    Output           : ${params.outdir}
    Annotations      : ${params.annotations_dir}
    Models           : ${params.models_dir}
    ──────────────────────────────────────────────────────────────────────────
    """.stripIndent()

    // Validate
    if (!params.input)           { error "Please provide --input" }
    if (!params.annotations_dir) { error "Please provide --annotations_dir" }
    if (!params.models_dir)      { error "Please provide --models_dir" }

    // Input channel
    if (params.input.endsWith('.csv')) {
        ch_sv_bed = Channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row -> [ [id: row.sample], file(row.bed, checkIfExists: true) ] }
    } else {
        ch_sv_bed = Channel
            .fromPath(params.input, checkIfExists: true)
            .map { bed -> [ [id: bed.baseName], bed ] }
    }

    // Reference channels
    ch_annotations = Channel.fromPath(params.annotations_dir, checkIfExists: true, type: 'dir').first()
    ch_models      = Channel.fromPath(params.models_dir, checkIfExists: true, type: 'dir').first()
    ch_scripts     = Channel.fromPath("${projectDir}/bin", checkIfExists: true, type: 'dir').first()
    ch_genome      = Channel.fromPath("${projectDir}/assets/caddsv/hg38.genome", checkIfExists: true).first()
    ch_header      = Channel.fromPath("${projectDir}/assets/caddsv/annotation_header.txt", checkIfExists: true).first()

    // Run pipeline
    CADDSV (
        ch_sv_bed,
        ch_annotations,
        ch_models,
        ch_scripts,
        ch_genome,
        ch_header
    )

    // Output
    CADDSV.out.scored.view { meta, file -> "Scored: ${meta.id} -> ${file}" }
}