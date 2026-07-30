// ============================================================
// 00_preprocess_reads.nf
// ONT read preprocessing: porechop adapter removal, QC, chopper split
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: PORECHOP
// Remove adapters from ONT reads
// ============================================================
process PORECHOP {
    label 'reads_env'
    tag   "porechop_${reads.simpleName}"
    cpus  { 8 * task.attempt }
    memory { 16.GB * task.attempt }
    time   { 6.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::porechop'
    container 'quay.io/biocontainers/porechop:0.2.4--py311he264feb_9'
    publishDir "${params.outdir}/00_preprocessing/porechop", mode: params.publish_mode

    input:
    path reads

    output:
    path "ONT_reads_pc.fastq.gz", emit: reads_pc

    script:
    """
    porechop \\
      --input ${reads} \\
      -o ONT_reads_pc.fastq.gz \\
      --discard_middle
    """
    // """ //implement_abi wf after final round of tests
    // porechop_abi \\
    //     -i ${reads} \\
    //     --threads ${task.cpus} \\
    //     -abi \\
    //     --discard_middle
    // """
    stub:
    """
    touch ONT_reads_pc.fastq.gz
    echo "PORECHOP stub completed"
    """
}

// ============================================================
// Process: READ_QC
// Quality control: FastQC, MultiQC, NanoPlot
// ============================================================
process FASTQC {
    label 'reads_env'
    tag   { "read_qc_${stage}" }
    cpus  { 8 * task.attempt }
    memory { 16.GB * task.attempt }
    time   { 6.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::fastqc'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir { "${params.outdir}/00_preprocessing/qc/${stage}/fastqc" }, mode: params.publish_mode

    input:
    path reads_pc
    val  stage

    output:
    path "qc_fastqc", emit: qc_dir

    script:
    """
    mkdir -p qc_fastqc
    fastqc -t ${task.cpus}   ${reads_pc} -o qc_fastqc
    """
    stub:
    """
    mkdir -p qc_fastqc
    touch qc_fastqc/${reads_pc}.testrun.zip
    touch qc_fastqc/${reads_pc}.testrun.html
    """
}

process NANOPLOT {
    label 'reads_env'
    tag   { "read_qc_${stage}" }
    cpus  { 8 * task.attempt }
    memory { 16.GB * task.attempt }
    time   { 6.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'nanoplot'
    container 'quay.io/biocontainers/nanoplot:1.47.1--pyhdfd78af_0'
    publishDir { "${params.outdir}/00_preprocessing/qc/${stage}/nanoplot" }, mode: params.publish_mode

    input:
    path reads_pc
    val  stage

    output:
    path "qc_nanoplot", emit: qc_dir

    script:
    """
    mkdir -p qc_nanoplot
    NanoPlot --fastq ${reads_pc} --plots dot --legacy hex --N50 -o qc_nanoplot -t ${task.cpus} 
    """
    stub:
    """
    mkdir -p qc_nanoplot
    touch qc_nanoplot/${reads_pc}.testrun.txt
    """
}

process NANOSTAT {
    label 'reads_env'
    tag   { "read_qc_${stage}" }
    cpus  { 8 * task.attempt }
    memory { 16.GB * task.attempt }
    time   { 6.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'nanostat'
    container 'quay.io/biocontainers/nanostat:0.1.5--py35_0'
    publishDir { "${params.outdir}/00_preprocessing/qc/${stage}/nanostat" }, mode: params.publish_mode

    input:
    path reads_pc
    val  stage

    output:
    path "qc_nanostat", emit: qc_dir

    script:
    """
    mkdir -p qc_nanostat
    NanoStat --fastq ${reads_pc} -t ${task.cpus}  > qc_nanostat/nanostat.txt
    """
    stub:
    """
    mkdir -p qc_nanostat
    touch qc_nanostat/${reads_pc}.testrun.txt
    """
}

process MULTIQC {
    label 'reads_env'
    tag   { "read_qc_${stage}" }
    cpus  { 8 * task.attempt }
    memory { 16.GB * task.attempt }
    time   { 6.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::multiqc'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir { "${params.outdir}/00_preprocessing/qc/${stage}/multiqc" }, mode: params.publish_mode

    input:
    path fastqc_dirs,   stageAs: 'fastqc_input_*'
    path nanoplot_dirs, stageAs: 'nanoplot_input_*'
    path nanostat_dirs, stageAs: 'nanostat_input_*'
    val  stage

    output:
    path "qc", emit: qc_dir

    script:
    """
    multiqc fastqc_input_* nanoplot_input_* nanostat_input_* -o qc --interactive #change to '.' ? 
    """
    stub:
    """
    mkdir -p qc
    touch qc/multiqc_report.${stage}.testrun.html
    """
}

// ============================================================
// Process: CHOPPER_SPLIT
// Split reads into assembly (longer) and polishing (shorter) sets
// Optionally filter symbiont contamination
// ============================================================
process CHOPPER_SPLIT {
    label 'reads_env'
    tag   "chopper_split"
    cpus  { 8 * task.attempt }
    memory { 32.GB * task.attempt }
    time   { 12.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::chopper'
    container 'quay.io/biocontainers/chopper:0.13.0--h7f49ad2_0'
    publishDir "${params.outdir}/00_preprocessing/reads", mode: params.publish_mode

    input:
    path reads_pc
    path symbiont_ref  // optional; pass [] when not provided

    output:
    path "assembly_reads.fastq.gz",  emit: assembly_reads
    path "polishing_reads.fastq.gz", emit: polishing_reads

    script:
    def contam_opt = params.skip_symbiont_filter || !symbiont_ref
        ? ''
        : "--contam ${symbiont_ref}"
    def half_cpus = Math.max(1, (task.cpus / 2) as int)
    """
    gunzip -c ${reads_pc} \\
        | chopper -q 3 -l 1000 ${contam_opt} --threads ${half_cpus} \\
        | gzip > assembly_reads.fastq.gz &
    PID1=\$!

    gunzip -c ${reads_pc} \\
        | chopper -q 5 -l 500 ${contam_opt} --threads ${half_cpus} \\
        | gzip > polishing_reads.fastq.gz &
    PID2=\$!

    wait \$PID1 \$PID2
    """

    stub:
    """
    touch assembly_reads.fastq.gz polishing_reads.fastq.gz
    echo "CHOPPER_SPLIT stub completed"
    """
}

// // ============================================================
// // Workflow: preprocess_reads
// // Orchestrates porechop → QC → chopper split
// // ============================================================
// workflow preprocess_reads {

//     take:
//         ont_raw_reads
//         symbiont_ref

//     main:
//         PORECHOP(ont_raw_reads)
//         READ_QC(PORECHOP.out.reads_pc)
//         CHOPPER_SPLIT(PORECHOP.out.reads_pc, symbiont_ref)

//     emit:
//         reads_pc        = PORECHOP.out.reads_pc
//         assembly_reads  = CHOPPER_SPLIT.out.assembly_reads
//         polishing_reads = CHOPPER_SPLIT.out.polishing_reads
//         qc_dir          = READ_QC.out.qc_dir
// }
