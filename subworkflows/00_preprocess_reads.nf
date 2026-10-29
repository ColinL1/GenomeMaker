// ============================================================
// 00_preprocess_reads.nf
// ONT read preprocessing: porechop adapter removal, QC, chopper split
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: PORECHOP
// Remove adapters from ONT reads
// ============================================================
include { PORECHOP; CHOPPER_SPLIT; CHOPPER_HQ_READS } from '../modules/00_preprocess_reads.nf'
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
        // Accumulate tool versions. Conditional stages mix INSIDE their branch so a
        // skipped stage contributes nothing rather than referencing a process that
        // never ran.
        def versions_ch = channel.empty()

        READ_QC(ont_raw_reads, 'raw')
        versions_ch = versions_ch.mix(READ_QC.out.versions)

        if (!params.run_porechop.toString().toBoolean()) {
            log.info "Skipping PORECHOP (run_porechop=false); using raw reads directly"
            trimmed_reads  = ont_raw_reads
            trimmed_qc_dir = READ_QC.out.qc_dir
            trimmed_raw_qc = READ_QC.out.fastqc_dir
                .mix(READ_QC.out.nanoplot_dir, READ_QC.out.nanostat_dir)
        } else {
            PORECHOP(ont_raw_reads)
            versions_ch = versions_ch.mix(PORECHOP.out.versions)
            trimmed_reads = PORECHOP.out.reads_pc
            READ_QC_TRIMMED(trimmed_reads, 'trimmed')
            trimmed_qc_dir = READ_QC_TRIMMED.out.qc_dir
            trimmed_raw_qc = READ_QC_TRIMMED.out.fastqc_dir
                .mix(READ_QC_TRIMMED.out.nanoplot_dir, READ_QC_TRIMMED.out.nanostat_dir)
            versions_ch = versions_ch.mix(READ_QC_TRIMMED.out.versions)
        }

        CHOPPER_SPLIT(trimmed_reads, symbiont_ref)

        // Separate high-quality subset for k-mer genome sizing -- see the
        // CHOPPER_HQ_READS comment for why this is not a CHOPPER_SPLIT output.
        CHOPPER_HQ_READS(trimmed_reads)

        versions_ch = versions_ch.mix(CHOPPER_SPLIT.out.versions, CHOPPER_HQ_READS.out.versions)

    emit:
        reads_pc         = trimmed_reads
        assembly_reads   = CHOPPER_SPLIT.out.assembly_reads
        polishing_reads  = CHOPPER_SPLIT.out.polishing_reads
        kmer_reads       = CHOPPER_HQ_READS.out.hq_reads
        versions         = versions_ch
        qc_dir_raw       = READ_QC.out.qc_dir
        qc_dir_trimmed   = trimmed_qc_dir
        // Raw QC directories for the final MultiQC recap (see 00b_qc_reads.nf).
        raw_qc_files     = READ_QC.out.fastqc_dir
            .mix(READ_QC.out.nanoplot_dir, READ_QC.out.nanostat_dir)
            .mix(trimmed_raw_qc)
}
