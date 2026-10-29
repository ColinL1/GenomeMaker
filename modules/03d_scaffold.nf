// ============================================================
// 03d_scaffold.nf
// Reference-guided scaffolding with RagTag, AFTER purge_dups and the contamination
// screens: haplotigs are purged while they are still separate contigs (RagTag would
// otherwise place both haplotypes side by side inside a chromosome), and
// contaminant contigs are never placed into chromosomes. Then the name map.
// ============================================================

nextflow.enable.dsl = 2

process RAGTAG_SCAFFOLD {
    tag "ragtag"

    conda 'bioconda::ragtag=2.1.0'
    container 'quay.io/biocontainers/ragtag:2.1.0--pyhdfd78af_2'
    publishDir "${params.outdir}/03d_scaffolding", mode: params.publish_mode

    input:
    path contigs
    path reference

    output:
    path "ragtag.scaffold.fasta", emit: scaffolds
    path "ragtag.scaffold.agp",   emit: agp
    path "ragtag.scaffold.stats", emit: stats
    path "ragtag.scaffold.confidence.txt"
    path "versions.yml", emit: versions

    script:
    // Defaults (minimap2 -x asm5, 100 bp unknown gaps) unless ragtag_args is set;
    // -r infers gap sizes from the reference.
    def extra = [params.ragtag_infer_gaps.toString().toBoolean() ? '-r' : '', params.ragtag_args ?: ''].findAll().join(' ')
    """
    ragtag.py scaffold -t ${task.cpus} ${extra} -o out ${reference} ${contigs} > ragtag.log 2>&1
    cp out/ragtag.scaffold.fasta out/ragtag.scaffold.agp out/ragtag.scaffold.stats out/ragtag.scaffold.confidence.txt .
    test -s ragtag.scaffold.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ragtag: "\$(ragtag.py --version 2>&1 | sed 's/^v//')"
        minimap2: "\$(minimap2 --version 2>&1)"
        reference: "${reference.name}"
    END_VERSIONS
    """

    stub:
    """
    cp ${contigs} ragtag.scaffold.fasta
    touch ragtag.scaffold.agp ragtag.scaffold.stats ragtag.scaffold.confidence.txt
    printf '"${task.process}":\\n    ragtag: stub\\n' > versions.yml
    """
}

// Final scaffold names -> RagTag scaffold names (by sequence MD5; FUNANNOTATE_SORT only renames):
// which scaffold_N is which reference-chromosome scaffold (<ref>_RagTag). Run only with RagTag.
process SCAFFOLD_NAME_MAP {
    label 'res_light'
    tag   "scaffold_name_map"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/03d_scaffolding/sorted", mode: params.publish_mode

    input:
    path presort, stageAs: 'presort/*'
    path final_fasta, stageAs: 'final/*'

    output:
    path "scaffold_name_map.tsv", emit: map

    script:
    """
    scaffold_name_map.py ${presort} ${final_fasta} scaffold_name_map.tsv
    """

    stub:
    """
    touch scaffold_name_map.tsv
    """
}
