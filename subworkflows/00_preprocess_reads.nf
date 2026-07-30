// ============================================================
// 00_preprocess_reads.nf
// ONT read preprocessing: porechop adapter removal, QC, chopper split
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: PORECHOP
// Remove adapters from ONT reads
// ============================================================
include { PORECHOP; CHOPPER_SPLIT } from '../modules/00_preprocess_reads.nf'
include { READ_QC; READ_QC as READ_QC_TRIMMED } from '../subworkflows/00b_qc_reads.nf'

// ============================================================
// Workflow: PREPROCESS_READS
// Orchestrates QC (raw) → porechop → QC (trimmed) → chopper split
// ============================================================
workflow PREPROCESS_READS {
    take:
        ont_raw_reads
        symbiont_ref

    main:
        READ_QC(ont_raw_reads, 'raw')

        if (params.skip_porechop) {
            log.info "Skipping PORECHOP (skip_porechop=true); using raw reads directly"
            trimmed_reads  = ont_raw_reads
            trimmed_qc_dir = READ_QC.out.qc_dir
        } else {
            PORECHOP(ont_raw_reads)
            trimmed_reads = PORECHOP.out.reads_pc
            READ_QC_TRIMMED(trimmed_reads, 'trimmed')
            trimmed_qc_dir = READ_QC_TRIMMED.out.qc_dir
        }

        CHOPPER_SPLIT(trimmed_reads, symbiont_ref)

    emit:
        reads_pc         = trimmed_reads
        assembly_reads   = CHOPPER_SPLIT.out.assembly_reads
        polishing_reads  = CHOPPER_SPLIT.out.polishing_reads
        qc_dir_raw       = READ_QC.out.qc_dir
        qc_dir_trimmed   = trimmed_qc_dir
}
