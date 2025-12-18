/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CADDSV SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CADDSV_PREPBED                          } from '../../../modules/local/caddsv/prepbed/main'
include { CADDSV_ANNOTATE as CADDSV_ANNOTATE_MAIN } from '../../../modules/local/caddsv/annotate/main'
include { CADDSV_ANNOTATE as CADDSV_ANNOTATE_UP   } from '../../../modules/local/caddsv/annotate/main'
include { CADDSV_ANNOTATE as CADDSV_ANNOTATE_DOWN } from '../../../modules/local/caddsv/annotate/main'
include { CADDSV_SCORE                            } from '../../../modules/local/caddsv/score/main'

workflow CADDSV {

    take:
    ch_sv_bed           // channel: [ val(meta), path(bed) ]
    ch_annotations_dir  // channel: path(annotations)
    ch_models_dir       // channel: path(models)
    ch_scripts_dir      // channel: path(scripts)
    ch_genome_file      // channel: path(genome)
    ch_header_file      // channel: path(header)

    main:

    ch_versions = Channel.empty()

    //
    // STEP 1: Prepare BED files (creates wchr/nochr/merged/flank versions)
    //
    CADDSV_PREPBED (
        ch_sv_bed,
        ch_genome_file
    )
    ch_versions = ch_versions.mix(CADDSV_PREPBED.out.versions)

    //
    // STEP 2A: Annotate main SV regions
    // ANNOTATE expects: tuple val(meta), path(bed_wchr), path(bed_nochr), path(merged_wchr), path(merged_nochr)
    //
    ch_main_for_annotate = CADDSV_PREPBED.out.bed_wchr
        .join(CADDSV_PREPBED.out.bed_nochr)
        .join(CADDSV_PREPBED.out.merged_wchr)
        .join(CADDSV_PREPBED.out.merged_nochr)

    CADDSV_ANNOTATE_MAIN (
        ch_main_for_annotate,
        ch_annotations_dir,
        ch_scripts_dir
    )
    ch_versions = ch_versions.mix(CADDSV_ANNOTATE_MAIN.out.versions.first())

    //
    // STEP 2B: Annotate 100bp upstream flanks
    // For flanks, we use the flank files for both regular and merged (flanks are small)
    //
    ch_up_for_annotate = CADDSV_PREPBED.out.flank_up_wchr
        .join(CADDSV_PREPBED.out.flank_up_nochr)
        .map { meta, wchr, nochr -> 
            // Use same files for merged since flanks are small non-overlapping regions
            [ meta, wchr, nochr, wchr, nochr ]
        }

    CADDSV_ANNOTATE_UP (
        ch_up_for_annotate,
        ch_annotations_dir,
        ch_scripts_dir
    )

    //
    // STEP 2C: Annotate 100bp downstream flanks
    //
    ch_down_for_annotate = CADDSV_PREPBED.out.flank_down_wchr
        .join(CADDSV_PREPBED.out.flank_down_nochr)
        .map { meta, wchr, nochr -> 
            [ meta, wchr, nochr, wchr, nochr ]
        }

    CADDSV_ANNOTATE_DOWN (
        ch_down_for_annotate,
        ch_annotations_dir,
        ch_scripts_dir
    )

    //
    // STEP 3: Score SVs - combine all three annotation sets
    //
    ch_for_scoring = CADDSV_ANNOTATE_MAIN.out.annotated
        .join(CADDSV_ANNOTATE_UP.out.annotated)
        .join(CADDSV_ANNOTATE_DOWN.out.annotated)

    CADDSV_SCORE (
        ch_for_scoring,
        ch_models_dir,
        ch_scripts_dir,
        ch_header_file
    )
    ch_versions = ch_versions.mix(CADDSV_SCORE.out.versions.first())

    emit:
    scored       = CADDSV_SCORE.out.scored
    scored_phred = CADDSV_SCORE.out.scored_phred
    annotated    = CADDSV_ANNOTATE_MAIN.out.annotated
    versions     = ch_versions
}