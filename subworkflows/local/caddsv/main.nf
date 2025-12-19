/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CADDSV SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CADDSV_PREPBED                          } from '../../../modules/local/caddsv/prepbed/main'
include { CADDSV_IDBED                            } from '../../../modules/local/caddsv/idbed/main'
include { CADDSV_ANNOTATE as CADDSV_ANNOTATE_MAIN } from '../../../modules/local/caddsv/annotate/main'
include { CADDSV_ANNOTATE as CADDSV_ANNOTATE_UP   } from '../../../modules/local/caddsv/annotate/main'
include { CADDSV_ANNOTATE as CADDSV_ANNOTATE_DOWN } from '../../../modules/local/caddsv/annotate/main'
include { CADDSV_ASSEMBLE                         } from '../../../modules/local/caddsv/assemble/main'
include { CADDSV_SCORE                            } from '../../../modules/local/caddsv/score/main'  // keep commented for now

workflow CADDSV {

    take:
    ch_sv_bed
    ch_annotations_dir
    ch_models_dir
    ch_scripts_dir
    ch_genome_file
    ch_header_file

    main:
    ch_versions = Channel.empty()

    CADDSV_IDBED(ch_sv_bed)
    ch_versions = ch_versions.mix(CADDSV_IDBED.out.versions)

    CADDSV_PREPBED(ch_sv_bed, ch_genome_file)
    ch_versions = ch_versions.mix(CADDSV_PREPBED.out.versions)

    CADDSV_ANNOTATE_MAIN(CADDSV_PREPBED.out.bed_main, ch_annotations_dir, ch_scripts_dir)
    ch_versions = ch_versions.mix(CADDSV_ANNOTATE_MAIN.out.versions.first())

    CADDSV_ANNOTATE_UP(CADDSV_PREPBED.out.bed_flank_up, ch_annotations_dir, ch_scripts_dir)
    CADDSV_ANNOTATE_DOWN(CADDSV_PREPBED.out.bed_flank_down, ch_annotations_dir, ch_scripts_dir)

    ch_annot_triplet = CADDSV_ANNOTATE_MAIN.out.annotated
        .join(CADDSV_ANNOTATE_UP.out.annotated)
        .join(CADDSV_ANNOTATE_DOWN.out.annotated)

    CADDSV_ASSEMBLE(
        ch_annot_triplet,
        ch_header_file
    )
    ch_versions = ch_versions.mix(CADDSV_ASSEMBLE.out.versions.first())

    // Join IDBED with assembled matrices on meta
    ch_score_input = CADDSV_IDBED.out.idbed.join(
        CADDSV_ASSEMBLE.out.matrix_main
            .join(CADDSV_ASSEMBLE.out.matrix_up)
            .join(CADDSV_ASSEMBLE.out.matrix_down)
    )

    // ch_score_input:
    // (meta, idbed) + (meta, matrix_main) + (meta, matrix_up) + (meta, matrix_down)
    ch_score_input = CADDSV_IDBED.out.idbed
        .join(CADDSV_ASSEMBLE.out.matrix_main)
        .join(CADDSV_ASSEMBLE.out.matrix_up)
        .join(CADDSV_ASSEMBLE.out.matrix_down)
        .map { meta, idbed, matrix_main, matrix_up, matrix_down ->
            tuple(meta, idbed, matrix_main, matrix_up, matrix_down)
        }

    CADDSV_SCORE(
        ch_score_input,
        ch_models_dir,
        ch_scripts_dir,
        ch_genome_file
    )
    
    ch_versions = ch_versions.mix(CADDSV_SCORE.out.versions.first())

    emit:
    idbed       = CADDSV_IDBED.out.idbed
    matrix_main = CADDSV_ASSEMBLE.out.matrix_main
    matrix_up   = CADDSV_ASSEMBLE.out.matrix_up
    matrix_down = CADDSV_ASSEMBLE.out.matrix_down
    annotated   = CADDSV_ANNOTATE_MAIN.out.annotated
    versions    = ch_versions
}