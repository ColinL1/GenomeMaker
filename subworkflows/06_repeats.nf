// ============================================================
// 06_repeats.nf
// Repeat identification, masking, and annotation
// barrnap rRNA → RepeatModeler → EDTA → combine → RepeatMasker
// ============================================================

nextflow.enable.dsl = 2

include { BARRNAP_RRNA ; REPEATMODELER ; EDTA_RUN ; COMBINE_REPEAT_DB ; REPEATMASKER_RUN } from '../modules/06_repeats.nf'
include { REPEAT_CLASSIFY ; REPEAT_GFF } from '../modules/06b_repeat_classify.nf'

// ============================================================
// Workflow: REPEATS_PIPELINE
// rRNA prediction → RepeatModeler + EDTA → combine → RepeatMasker
// ============================================================
workflow REPEATS_PIPELINE {

    take:
        assembly

    main:
        BARRNAP_RRNA(assembly)
        REPEATMODELER(assembly)
        EDTA_RUN(assembly)
        COMBINE_REPEAT_DB(REPEATMODELER.out.repeat_lib, EDTA_RUN.out.edta_lib)
        // Unknown families classified BEFORE masking, so RepeatMasker reports final labels
        REPEAT_CLASSIFY(COMBINE_REPEAT_DB.out.dedup_lib, params.rm_libraries ?: '')
        REPEATMASKER_RUN(assembly, REPEAT_CLASSIFY.out.lib)
        REPEAT_GFF(REPEATMASKER_RUN.out.repeat_tbl)

    emit:
        masked_assembly = REPEATMASKER_RUN.out.masked_assembly
        repeat_tbl      = REPEATMASKER_RUN.out.repeat_tbl
        repeat_summary  = REPEATMASKER_RUN.out.repeat_summary
        repeat_lib      = REPEAT_CLASSIFY.out.lib           // classified library used for masking
        repeat_gff      = REPEAT_GFF.out.gff
        classify_report = REPEAT_CLASSIFY.out.report
        rrna_gff        = BARRNAP_RRNA.out.rrna_gff
        versions        = BARRNAP_RRNA.out.versions.mix(REPEATMODELER.out.versions, EDTA_RUN.out.versions, COMBINE_REPEAT_DB.out.versions, REPEAT_CLASSIFY.out.versions, REPEATMASKER_RUN.out.versions)
}
