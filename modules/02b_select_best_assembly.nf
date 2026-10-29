// ============================================================
// 02b_select_best_assembly.nf
// Rank assemblers (bin/rank.py): size outliers last; within the BUSCO tie band
// the assembly closest to the expected size wins; N50 breaks a size tie.
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
    path quast_reports       // ${name}_quast_report/ dirs
    val  expected_genome_size

    output:
    path "best_assembler.txt",   emit: best_name
    path "assembly_ranking.tsv", emit: ranking_table
    path "versions.yml", emit: versions

    script:
    // Expected haploid size (params.genome_size or GenomeScope2); required by rank.py.
    def expected = (expected_genome_size ?: 0)
    """
    export EXPECTED_GENOME_SIZE='${expected}'
    export SIZE_TOLERANCE='${params.assembly_size_tolerance}'
    export BUSCO_TIE_BAND='${params.busco_tie_band}'
    rank.py
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "\$(python3 --version 2>&1 | sed 's/Python //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    echo "flye" > best_assembler.txt
    printf 'assembler\tcomplete_pct\tduplicated_pct\tN50\ttotal_bp\tsize_dev\tsize_check\tquast_contigs\tquast_largest_contig\tselected\tdecided_by\n' > assembly_ranking.tsv
    printf 'flye\t98.0\t1.0\t1000000\t600000000\t+0.000\tok\t\t\tyes\tstub\n' >> assembly_ranking.tsv
    echo "RANK_ASSEMBLIES stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
