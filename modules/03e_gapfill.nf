// ============================================================
// 03e_gapfill.nf
// Gap filling of the RagTag scaffolds (params.run_gapfill): TGS-GapCloser --ne with the
// porechopped ONT reads, then a second Medaka pass (MEDAKA_POLISH as MEDAKA_GAPFILL, main.nf).
// No gate: the gap-filled (and polished) scaffolds go on. The AGP is rebuilt from the
// N runs left, and the report counts what each step closed.
// ============================================================

nextflow.enable.dsl = 2

// TGS-GapCloser reads FASTA only
process READS_TO_FASTA {
    label 'res_light'
    tag   "reads_to_fasta"

    conda 'bioconda::seqkit'
    container 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_0'

    input:
    path reads

    output:
    path "reads.fa",     emit: fasta
    path "versions.yml", emit: versions

    script:
    """
    gzip -dc ${reads} | seqkit fq2fa -o reads.fa
    test -s reads.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: "\$(seqkit version 2>&1 | sed 's/seqkit v//' || echo unknown)"
    END_VERSIONS
    """

    stub:
    """
    touch reads.fa
    printf '"${task.process}":\\n    seqkit: stub\\n' > versions.yml
    """
}

// --ne: no error correction of the candidate reads (no racon); fills come from the raw reads.
// Only the read mapping is multi-threaded; the fill step is single-threaded and memory-bound.
// scratch: the .paf files (44 GB on the Porites assembly) stay on node-local disk.
process TGS_GAPFILL {
    tag "tgsgapcloser"
    scratch true

    conda 'bioconda::tgsgapcloser=1.2.1 bioconda::racon conda-forge::util-linux'
    // The biocontainer has no getopt (util-linux), which tgsgapcloser needs
    container 'oras://community.wave.seqera.io/library/tgsgapcloser_racon_util-linux:236f7f9102d6b420'
    publishDir "${params.outdir}/03e_gapfill/tgsgapcloser", mode: params.publish_mode, pattern: 'tgs.*'

    input:
    path scaffolds
    path reads

    output:
    path "tgs.scaff_seqs",      emit: scaffolds
    path "tgs.gap_fill_detail", emit: detail
    path "versions.yml",        emit: versions

    script:
    """
    tgsgapcloser --scaff ${scaffolds} --reads ${reads} --output \$PWD/tgs --ne --tgstype ont \\
      --min_nread ${params.gapfill_min_nread} --thread ${task.cpus} > tgs.log 2>&1
    test -s tgs.scaff_seqs
    rm -f tgs.*.paf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tgsgapcloser: "\$(tgsgapcloser 2>&1 | sed -n 's/.*Version : \\([^ ;]*\\).*/\\1/p' | head -1)"
        minimap2: "\$(minimap2 --version 2>&1)"
    END_VERSIONS
    """

    stub:
    """
    cp ${scaffolds} tgs.scaff_seqs
    touch tgs.gap_fill_detail
    printf '"${task.process}":\\n    tgsgapcloser: stub\\n' > versions.yml
    """
}

// AGP of the gap-filled scaffolds: one gap line per N run left, RagTag's gap type and linkage columns
process GAPFILL_AGP {
    label 'res_light'
    tag   "gapfill_agp"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/03e_gapfill", mode: params.publish_mode

    input:
    path scaffolds
    path ragtag_agp

    output:
    path "gapfilled.agp", emit: agp
    path "versions.yml",  emit: versions

    script:
    """
    agp_from_gaps.py ${scaffolds} ${ragtag_agp} gapfilled.agp

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version | sed 's/Python //')" > versions.yml
    """

    stub:
    """
    touch gapfilled.agp
    printf '"${task.process}":\\n    python: stub\\n' > versions.yml
    """
}

// MultiQC table: gaps and contig N50 after RagTag, TGS-GapCloser and Medaka; joins and fills
process GAPFILL_REPORT {
    label 'res_light'
    tag   "gapfill_report"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/03e_gapfill/qc", mode: params.publish_mode

    input:
    path detail
    path ragtag_agp
    path gfastats_ragtag, stageAs: 'g/ragtag/*'
    path gfastats_tgs,    stageAs: 'g/tgs/*'
    path gfastats_medaka, stageAs: 'g/medaka/*'   // or a NO_ placeholder (run_gapfill_polish=false)

    output:
    path "gapfill_mqc.tsv", emit: mqc
    path "versions.yml",    emit: versions

    script:
    """
    gapfill_report.py ${detail} ${ragtag_agp} ${gfastats_ragtag} ${gfastats_tgs} ${gfastats_medaka} gapfill_mqc.tsv

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version | sed 's/Python //')" > versions.yml
    """

    stub:
    """
    printf "# id: 'gapfill'\\n# plot_type: 'table'\\nStage\\tGaps\\n" > gapfill_mqc.tsv
    printf '"${task.process}":\\n    python: stub\\n' > versions.yml
    """
}
