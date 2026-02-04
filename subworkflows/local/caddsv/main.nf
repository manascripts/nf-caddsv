/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CADDSV SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CADDSV_PREPARE_REFERENCES               } from '../../../modules/local/caddsv/prepare_references/main'
include { CADDSV_PREPBED                          } from '../../../modules/local/caddsv/prepbed/main'
include { CADDSV_ANNOTATE as CADDSV_ANNOTATE_MAIN } from '../../../modules/local/caddsv/annotate/main'
include { CADDSV_ANNOTATE as CADDSV_ANNOTATE_UP   } from '../../../modules/local/caddsv/annotate/main'
include { CADDSV_ANNOTATE as CADDSV_ANNOTATE_DOWN } from '../../../modules/local/caddsv/annotate/main'
include { CADDSV_SCORE                            } from '../../../modules/local/caddsv/score/main'

workflow CADDSV {

    take:
    ch_sv_bed
    ch_annotations_dir
    ch_models_dir
    ch_scripts_dir
    ch_genome_file

    main:
    ch_versions = Channel.empty()

    if (!params.skip_prepare_references) {
        CADDSV_PREPARE_REFERENCES()
        ch_versions = ch_versions.mix(CADDSV_PREPARE_REFERENCES.out.versions)
    }

    CADDSV_PREPBED(ch_sv_bed, ch_genome_file)
    ch_versions = ch_versions.mix(CADDSV_PREPBED.out.versions)

    CADDSV_ANNOTATE_MAIN(CADDSV_PREPBED.out.bed_main, ch_annotations_dir, ch_scripts_dir)
    ch_versions = ch_versions.mix(CADDSV_ANNOTATE_MAIN.out.versions.first())

    CADDSV_ANNOTATE_UP(CADDSV_PREPBED.out.bed_flank_up, ch_annotations_dir, ch_scripts_dir)
    CADDSV_ANNOTATE_DOWN(CADDSV_PREPBED.out.bed_flank_down, ch_annotations_dir, ch_scripts_dir)

    CADDSV_ANNOTATE_MAIN.out.annotated
        .join(CADDSV_ANNOTATE_UP.out.annotated)
        .join(CADDSV_ANNOTATE_DOWN.out.annotated)

    // Build one tuple per sample: (meta, original_bed, matrix_main, matrix_up, matrix_down)
    ch_score_input = CADDSV_PREPBED.out.sorted_original
        .join(CADDSV_ANNOTATE_MAIN.out.annotated)
        .join(CADDSV_ANNOTATE_UP.out.annotated)
        .join(CADDSV_ANNOTATE_DOWN.out.annotated)
        .map { meta, original_bed, matrix_main, matrix_up, matrix_down ->
            tuple(meta, original_bed, matrix_main, matrix_up, matrix_down)
        }

    CADDSV_SCORE(
        ch_score_input,
        ch_models_dir,
        ch_scripts_dir,
        ch_genome_file, 
        ch_annotations_dir
    )
    ch_versions = ch_versions.mix(CADDSV_SCORE.out.versions.first())

    emit:
    score       = CADDSV_SCORE.out.score
    annotated   = CADDSV_ANNOTATE_MAIN.out.annotated
    versions    = ch_versions
}