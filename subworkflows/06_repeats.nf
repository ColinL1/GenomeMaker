// ============================================================
// 06_repeats.nf
// Repeat identification, masking, and annotation
// barrnap rRNA → RepeatModeler → EDTA → combine → RepeatMasker
// ============================================================

nextflow.enable.dsl = 2

include { BARRNAP_RRNA ; REPEATMODELER ; EDTA_RUN ; COMBINE_REPEAT_DB ; REPEATMASKER_RUN } from '../modules/06_repeats.nf'

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
        REPEATMASKER_RUN(assembly, COMBINE_REPEAT_DB.out.dedup_lib)

    emit:
        masked_assembly = REPEATMASKER_RUN.out.masked_assembly
        repeat_tbl      = REPEATMASKER_RUN.out.repeat_tbl
}
