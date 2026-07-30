// ============================================================
// 05_mitogenome.nf
// Mitochondrial genome extraction by mapping assembly reads
// to a reference mitogenome
// ============================================================

nextflow.enable.dsl = 2

include { MITO_READ_FILTER ; MITO_MAPPING ; MITO_COVERAGE_STATS } from '../modules/05_mitogenome.nf'

// ============================================================
// Workflow: MITOGENOME_PIPELINE
// Reads filter → map to reference → coverage stats
// ============================================================
workflow MITOGENOME_PIPELINE {

    take:
        assembly_reads
        mito_ref

    main:
        MITO_READ_FILTER(assembly_reads)
        MITO_MAPPING(MITO_READ_FILTER.out.mito_reads, mito_ref)
        MITO_COVERAGE_STATS(MITO_MAPPING.out.mito_bam, MITO_MAPPING.out.mito_bai)

    emit:
        mito_bam    = MITO_MAPPING.out.mito_bam
        mito_bai    = MITO_MAPPING.out.mito_bai
        coverage    = MITO_COVERAGE_STATS.out.coverage
        flagstat    = MITO_COVERAGE_STATS.out.flagstat
}
