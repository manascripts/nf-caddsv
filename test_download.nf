#!/usr/bin/env nextflow
/*
 * Simple test workflow - download CADD-SV references only
 */

nextflow.enable.dsl = 2

include { CADDSV_PREPARE_REFERENCES } from './modules/local/caddsv/prepare_references/main'

workflow {
    log.info """
    ══════════════════════════════════════════════════════════════════════════
      Testing CADD-SV Reference Download
    ══════════════════════════════════════════════════════════════════════════
    Cache dir: ${params.caddsv_cache_dir}
    ──────────────────────────────────────────────────────────────────────────
    """.stripIndent()

    CADDSV_PREPARE_REFERENCES()

    CADDSV_PREPARE_REFERENCES.out.annotations_dir.view { "Annotations: $it" }
    CADDSV_PREPARE_REFERENCES.out.models_dir.view { "Models: $it" }
}