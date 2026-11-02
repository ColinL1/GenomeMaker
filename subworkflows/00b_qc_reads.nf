// ============================================================
// 00b_qc_reads.nf
// Read QC: FastQC, NanoPlot, NanoStat, MultiQC
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
        fastqc_dir   = FASTQC.out.qc_dir
        nanoplot_dir = NANOPLOT.out.qc_dir
        nanostat_dir = NANOSTAT.out.qc_dir
        versions        = FASTQC.out.versions.mix(NANOPLOT.out.versions, NANOSTAT.out.versions, MULTIQC.out.versions)
}
