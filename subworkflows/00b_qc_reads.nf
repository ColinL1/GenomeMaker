// ============================================================
// 00_preprocess_reads.nf
// ONT read preprocessing: porechop adapter removal, QC, chopper split
// ============================================================

nextflow.enable.dsl = 2

include { FASTQC; NANOPLOT; NANOSTAT; MULTIQC } from '../modules/00_preprocess_reads.nf'

// ============================================================
// Workflow: READ_QC
// ============================================================
workflow READ_QC {
    take:
        ont_raw_reads
        stage
    main:
        FASTQC(ont_raw_reads, stage)
        NANOPLOT(ont_raw_reads, stage)
        NANOSTAT(ont_raw_reads, stage)
        MULTIQC(
            FASTQC.out.qc_dir.collect(),
            NANOPLOT.out.qc_dir.collect(),
            NANOSTAT.out.qc_dir.collect(),
            stage
        )
    emit:
        qc_dir       = MULTIQC.out.qc_dir
}
