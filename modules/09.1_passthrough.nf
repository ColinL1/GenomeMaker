// ============================================================
// 99_passthrough.nf
// Pass-through process: copy input file to output unchanged
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: COPY_FILE
// Copy a file to output (used when optional stages are skipped)
// ============================================================
process COPY_FILE {
    tag   "passthrough_${file.simpleName}"
    cpus  { 1 * task.attempt }
    memory { 2.GB * task.attempt }
    time   { 1.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    // conda 'bioconda::samtools'
    // container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'
    publishDir "${params.outdir}/06_repeats", mode: params.publish_mode

    input:
    path file

    output:
    path "${file.simpleName}", emit: output

    script:
    """
    cp ${file} ${file.simpleName}
    """

    stub:
    """
    touch ${file.simpleName}
    echo "COPY_FILE stub completed"
    """
}
