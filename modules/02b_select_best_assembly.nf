// ============================================================
// 02b_select_best_assembly.nf
// Rank assemblers by BUSCO completeness (penalizing duplication)
// with N50 as a tie-breaker, and select the winner.
// ============================================================

nextflow.enable.dsl = 2

process RANK_ASSEMBLIES {
    tag "rank_assemblies"
    cpus  { 2 * task.attempt }
    memory { 4.GB * task.attempt }
    time   { 30.min * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'conda-forge::python'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/02_assembly/selection", mode: params.publish_mode

    input:
    path gfastats_reports    // ${name}.gfastats.txt  (name-embedded filenames)
    path busco_summaries     // ${name}.${lineage}.short_summary.txt, HOST lineage only

    output:
    path "best_assembler.txt", emit: best_name
    path "assembly_ranking.tsv", emit: ranking_table

    script:
    """
    python3 ${params.bin_path}rank.py
    """

    stub:
    """
    echo "flye" > best_assembler.txt
    touch assembly_ranking.tsv
    echo "RANK_ASSEMBLIES stub completed"
    """
}
