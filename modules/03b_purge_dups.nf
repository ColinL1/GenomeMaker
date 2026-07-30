// ============================================================
// 03b_purge_dups.nf
// Standard purge_dups pipeline (read-depth + self-alignment based
// haplotig purging), using ONT coverage. Runs on the SELECTED
// assembly, before BlobToolKit filtering.
// ============================================================

nextflow.enable.dsl = 2

process PURGE_DUPS_COVERAGE {
    label 'purge_dups_env'
    tag   "purge_dups_coverage"
    cpus  { params.threads_default * task.attempt }
    memory { 32.GB * task.attempt }
    time   { 8.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::purge_dups bioconda::minimap2'
    container 'quay.io/biocontainers/purge_dups:1.2.6--h577a1d6_3'
    publishDir "${params.outdir}/03b_purge_dups", mode: params.publish_mode

    input:
    path assembly
    path reads   // ONT reads (assembly_reads)

    output:
    path "PB.base.cov", emit: base_cov
    path "PB.stat",     emit: stat
    path "cutoffs",     emit: cutoffs

    script:
    """
    minimap2 -x map-ont -t ${task.cpus} ${assembly} ${reads} | gzip -c > reads_vs_asm.paf.gz
    pbcstat reads_vs_asm.paf.gz
    calcuts PB.stat > cutoffs
    """

    stub:
    """
    touch PB.base.cov PB.stat cutoffs
    echo "PURGE_DUPS_COVERAGE stub completed"
    """
}

process PURGE_DUPS_RUN {
    label 'purge_dups_env'
    tag   "purge_dups_run"
    cpus  { params.threads_default * task.attempt }
    memory { 32.GB * task.attempt }
    time   { 8.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::purge_dups bioconda::minimap2'
    container 'quay.io/biocontainers/purge_dups:1.2.6--h577a1d6_3'
    publishDir "${params.outdir}/03b_purge_dups", mode: params.publish_mode

    input:
    path assembly
    path base_cov
    path cutoffs

    output:
    path "purged.fasta",    emit: purged_assembly
    path "haplotigs.fasta", emit: haplotigs
    path "dups.bed",        emit: dups_bed

    script:
    """
    split_fa ${assembly} > split.fa
    minimap2 -x asm5 -DP -t ${task.cpus} split.fa split.fa | gzip -c > split_self.paf.gz

    purge_dups -2 -T ${cutoffs} -c ${base_cov} split_self.paf.gz > dups.bed

    get_seqs -e dups.bed ${assembly}
    # get_seqs writes purged.fa and hap.fa by default
    mv purged.fa purged.fasta
    mv hap.fa haplotigs.fasta
    """

    stub:
    """
    touch purged.fasta haplotigs.fasta dups.bed
    echo "PURGE_DUPS_RUN stub completed"
    """
}
