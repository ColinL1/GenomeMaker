// ============================================================
// 05_mitogenome.nf
// Mitochondrial genome extraction by mapping assembly reads
// to a reference mitogenome
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: MITO_READ_FILTER
// Filter assembly reads for high quality (Q>=10) for mitogenome
// ============================================================
process MITO_READ_FILTER {
    label 'mito_env'
    tag   "mito_read_filter"
    cpus  { 8 * task.attempt }
    memory { 16.GB * task.attempt }
    time   { 6.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::chopper'
    container 'quay.io/biocontainers/chopper:0.13.0--h7f49ad2_0'
    publishDir "${params.outdir}/05_mitogenome", mode: params.publish_mode

    input:
    path reads

    output:
    path "mito_reads_q10.fastq.gz", emit: mito_reads

    script:
    """
    # Filter for quality >= 10 for mitogenome mapping
    gunzip -c ${reads} | \\
      chopper -q 10 | \\
      gzip > mito_reads_q10.fastq.gz
    """

    stub:
    """
    touch mito_reads_q10.fastq.gz
    echo "MITO_READ_FILTER stub completed"
    """
}

// ============================================================
// Process: MITO_MAPPING
// Map filtered reads to mitogenome reference, produce BAM
// ============================================================
process MITO_MAPPING {
    label 'mito_env'
    tag   "mito_mapping"
    cpus  { params.threads_default * task.attempt }
    memory { 32.GB * task.attempt }
    time   { 12.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::minimap2 samtools'
    container 'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:836eb07132d2de763323bbd4f1083b3fdf328759-0'
    publishDir "${params.outdir}/05_mitogenome", mode: params.publish_mode

    input:
    path reads
    path mito_ref

    output:
    path "mito_mapped.sorted.bam",  emit: mito_bam
    path "mito_mapped.sorted.bam.bai", emit: mito_bai

    script:
    """
    # Create minimap2 index for reference
    minimap2 -d mito_ref.mmi ${mito_ref}

    # Map reads to mitogenome reference
    minimap2 -ax map-ont mito_ref.mmi ${reads} | \\
      samtools view -bS - | \\
      samtools sort -@ ${task.cpus} -o mito_mapped.sorted.bam

    samtools index mito_mapped.sorted.bam
    """

    stub:
    """
    touch mito_mapped.sorted.bam mito_mapped.sorted.bam.bai
    echo "MITO_MAPPING stub completed"
    """
}

// ============================================================
// Process: MITO_COVERAGE_STATS
// Generate coverage statistics for mitogenome
// ============================================================
process MITO_COVERAGE_STATS {
    label 'mito_env'
    tag   "mito_coverage"
    cpus  { 2 * task.attempt }
    memory { 4.GB * task.attempt }
    time   { 1.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::samtools'
    container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'
    publishDir "${params.outdir}/05_mitogenome", mode: params.publish_mode

    input:
    path bam
    path bai

    output:
    path "mito_coverage.txt", emit: coverage
    path "mito_flagstat.txt", emit: flagstat

    script:
    """
    samtools coverage ${bam} > mito_coverage.txt
    samtools flagstat ${bam} > mito_flagstat.txt
    """

    stub:
    """
    touch mito_coverage.txt mito_flagstat.txt
    echo "MITO_COVERAGE_STATS stub completed"
    """
}

// // ============================================================
// // Workflow: mitogenome_pipeline
// // Reads filter → map to reference → coverage stats
// // ============================================================
// workflow mitogenome_pipeline {

//     take:
//         assembly_reads
//         mito_ref

//     main:
//         MITO_READ_FILTER(assembly_reads)
//         MITO_MAPPING(MITO_READ_FILTER.out.mito_reads, mito_ref)
//         MITO_COVERAGE_STATS(MITO_MAPPING.out.mito_bam, MITO_MAPPING.out.mito_bai)

//     emit:
//         mito_bam    = MITO_MAPPING.out.mito_bam
//         mito_bai    = MITO_MAPPING.out.mito_bai
//         coverage    = MITO_COVERAGE_STATS.out.coverage
//         flagstat    = MITO_COVERAGE_STATS.out.flagstat
// }
