// ============================================================
// 02b_select_best_assembly.nf
// Rank assemblers by BUSCO completeness (penalizing duplication)
// with N50 as a tie-breaker, and select the winner.
// ============================================================

nextflow.enable.dsl = 2

process RANK_ASSEMBLIES {
    tag "rank_assemblies"

    conda 'conda-forge::python'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/02_assembly/selection", mode: params.publish_mode

    input:
    path gfastats_reports    // ${name}.gfastats.txt  (name-embedded filenames)
    path busco_summaries     // ${name}.${lineage}.short_summary.txt, HOST lineage only
    path quast_reports       // ${name}_quast_report/ dirs -- QUAST was previously computed and then unused
    path rank_script
    val  expected_genome_size

    output:
    path "best_assembler.txt",   emit: best_name
    path "assembly_ranking.tsv", emit: ranking_table
    path "versions.yml", emit: versions

    script:
    // EXPECTED_GENOME_SIZE lets the ranking reject an assembly whose total length
    // is wildly off the GenomeScope2 estimate -- a doubled assembly can score well
    // on BUSCO completeness while being badly wrong. Empty/0 disables the check.
    def expected = (expected_genome_size ?: 0)
    """
    PYTHON=\$(command -v python3 || command -v python) || { echo "ERROR: no python interpreter found" >&2; exit 1; }
    export EXPECTED_GENOME_SIZE='${expected}'
    export SIZE_TOLERANCE='${params.assembly_size_tolerance}'
    export BUSCO_TIE_BAND='${params.busco_tie_band}'
    \$PYTHON ${rank_script}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "\$(python3 --version 2>&1 | sed 's/Python //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    echo "flye" > best_assembler.txt
    touch assembly_ranking.tsv
    echo "RANK_ASSEMBLIES stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
