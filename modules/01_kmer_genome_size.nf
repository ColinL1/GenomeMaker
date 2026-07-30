// ============================================================
// 01_kmer_genome_size.nf
// K-mer counting, histogram, and GenomeScope 2.0 estimation
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: MERYL_COUNT
// Count k-mers from reads
// ============================================================
process MERYL_COUNT {
    label 'kmer_env'
    tag   "meryl_count"
    cpus  { 8 * task.attempt }
    memory { 64.GB * task.attempt }
    time   { 24.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::meryl'
    container 'quay.io/biocontainers/meryl:1.2--he1b5a44_0'
    publishDir "${params.outdir}/01_kmer_genome_size/meryl", mode: params.publish_mode

    input:
    path reads_pc

    output:
    path "kmer_db.meryl", emit: meryl_db

    script:
    """
    meryl count k=${params.kmer_size} output kmer_db.meryl ${reads_pc}
    """

    stub:
    """
    mkdir -p kmer_db.meryl
    echo "MERYL_COUNT stub completed"
    """
}

// ============================================================
// Process: MERYL_HIST
// Generate k-mer frequency histogram
// ============================================================
process MERYL_HIST {
    label 'kmer_env'
    tag   "meryl_histogram"
    cpus  { 4 * task.attempt }
    memory { 16.GB * task.attempt }
    time   { 6.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::meryl'
    container 'quay.io/biocontainers/meryl:1.2--he1b5a44_0'
    publishDir "${params.outdir}/01_kmer_genome_size/meryl", mode: params.publish_mode

    input:
    path meryl_db

    output:
    path "kmer_db.meryl.hist", emit: meryl_hist

    script:
    """
    echo "Running MERYL_HIST"
    meryl histogram kmer_db.meryl > kmer_db.meryl.hist
    """

    stub:
    """
    touch kmer_db.meryl.hist
    echo "MERYL_HIST stub completed"
    """
}

// ============================================================
// Process: GENOMESCOPE
// Run GenomeScope 2.0 to estimate genome properties
// ============================================================
process GENOMESCOPE {
    label 'kmer_env'
    tag   "genomescope"
    cpus  { 2 * task.attempt }
    memory { 8.GB * task.attempt }
    time   { 2.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'conda-forge::r-base conda-forge::r-ggplot2'
    container 'quay.io/biocontainers/mulled-v2-5ff8b00c2d7f6173e034c115dfe295627ff99689:beb0ad5f49ec2904f79edffacd41bba38492e881-0'
    publishDir "${params.outdir}/01_kmer_genome_size/genomescope", mode: params.publish_mode

    input:
    path meryl_hist
    path genomescope_script

    output:
    path "genomescope_output_${params.kmer_size}", emit: genomescope_out

    script:
    """
    awk 'BEGIN{OFS=" "}{print \$1, \$2}' ${meryl_hist} > kmer_db.meryl.hist.clean

    Rscript ${genomescope_script} kmer_db.meryl.hist.clean ${params.kmer_size} 150 genomescope_output_${params.kmer_size}
    """
    stub:
    """
    mkdir -p genomescope_output_${params.kmer_size}
    touch genomescope_output_${params.kmer_size}/summary.txt
    echo "GENOMESCOPE stub completed"
    """
}

process PARSE_GENOME_SIZE {
    tag "parse_genomescope"
    label 'kmer_env'

    input:
    path genomescope_out

    output:
    env 'GENOME_SIZE', emit: genome_size

    script:
    """
    SUMMARY=\$(ls ${genomescope_out}/*summary.txt | head -n1)

    GENOME_SIZE=\$(grep -m1 'Genome Haploid Length' "\$SUMMARY" \\
        | sed -E 's/^Genome Haploid Length[[:space:]]+([0-9,]+) bp.*/\\1/' \\
        | tr -d ',')

    if [ -z "\$GENOME_SIZE" ]; then
        echo "ERROR: could not parse Genome Haploid Length from \$SUMMARY" >&2
        exit 1
    fi
    echo "Parsed genome size: \$GENOME_SIZE bp"
    """

    stub:
    """
    GENOME_SIZE=600000000
    echo "PARSE_GENOME_SIZE stub completed"
    """
}

// // ============================================================
// // Workflow: kmer_genome_size
// // Orchestrates meryl count → histogram → genomescope
// // Emits estimated genome size for use by assemblers
// // ============================================================
// workflow kmer_genome_size {

//     take:
//         reads_pc
//         genomescope_script

//     main:
//         MERYL_COUNT(reads_pc)
//         MERYL_HIST(MERYL_COUNT.out.meryl_db)
//         GENOMESCOPE(MERYL_HIST.out.meryl_hist, genomescope_script)

//     emit:
//         meryl_db        = MERYL_COUNT.out.meryl_db
//         meryl_hist      = MERYL_HIST.out.meryl_hist
//         genomescope_out = GENOMESCOPE.out.genomescope_out
// }
