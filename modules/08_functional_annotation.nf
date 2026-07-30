// ============================================================
// 08_functional_annotation.nf
// Functional annotation pipeline:
// Phobius → InterProScan (funannotate iprscan) → eggNOG-mapper → funannotate annotate
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: PHOBIUS_RUN
// Predict signal peptides and transmembrane helices
// ============================================================
process PHOBIUS_RUN {
    label 'functional_env'
    tag   "phobius"
    cpus  { 4 * task.attempt }
    memory { 8.GB * task.attempt }
    time   { 4.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::phobius'
    container 'interpro/phobius:1.01'
    publishDir "${params.outdir}/08_functional_annotation/phobius", mode: params.publish_mode

    input:
    path proteins

    output:
    path "phobius_results.txt", emit: phobius_results

    script:
    """
    phobius.pl ${proteins} > phobius_results.txt
    """

    stub:
    """
    touch phobius_results.txt
    echo "PHOBIUS_RUN stub completed"
    """
}

// ============================================================
// Process: IPRSCAN_RUN
// InterProScan via funannotate wrapper (Docker)
// ============================================================
process IPRSCAN_RUN {
    label 'functional_env'
    tag   "iprscan"
    cpus  { params.threads_funannotate * task.attempt }
    memory { 64.GB * task.attempt }
    time   { 48.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::funannotate interproscan-db'
    container 'quay.io/biocontainers/funannotate:1.8.17--pyhdfd78af_5'
    publishDir "${params.outdir}/08_functional_annotation/interproscan", mode: params.publish_mode

    input:
    path proteins

    output:
    path "iprscan.xml", emit: iprscan_xml

    script:
    """
    funannotate iprscan \\
      -i ${proteins} \\
      -m docker \\
      -c ${task.cpus} \\
      -o iprscan.xml
    """

    stub:
    """
    touch iprscan.xml
    echo "IPRSCAN_RUN stub completed"
    """
}

// ============================================================
// Process: EGGNOG_MAPPER
// Run eggNOG-mapper for orthology and GO annotations
// ============================================================
process EGGNOG_MAPPER {
    label 'functional_env'
    tag   "eggnog"
    cpus  { params.threads_funannotate * task.attempt }
    memory { 64.GB * task.attempt }
    time   { 24.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::eggnog-mapper'
    container 'quay.io/biocontainers/eggnog-mapper:2.1.15--pyhdfd78af_0'
    publishDir "${params.outdir}/08_functional_annotation/eggnog", mode: params.publish_mode

    input:
    path proteins
    val  funannotate_db

    output:
    path "eggnog.emapper.annotations", emit: eggnog_anno
    path "eggnog.emapper.seed_orthologs", emit: eggnog_ortho

    script:
    """
    emapper.py \\
      --cpu ${task.cpus} \\
      -m mmseqs \\
      --data_dir ${funannotate_db} \\
      -i ${proteins} \\
      -o eggnog \\
      --no_egeny
    """

    stub:
    """
    touch eggnog.emapper.annotations eggnog.emapper.seed_orthologs
    echo "EGGNOG_MAPPER stub completed"
    """
}

// ============================================================
// Process: FUNANNOTATE_ANNOTATE
// Final functional annotation synthesis
// ============================================================
process FUNANNOTATE_ANNOTATE {
    label 'functional_env'
    tag   "funannotate_annotate"
    cpus  { params.threads_funannotate * task.attempt }
    memory { 128.GB * task.attempt }
    time   { 48.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::funannotate'
    container 'quay.io/biocontainers/funannotate:1.8.17--pyhdfd78af_5'
    publishDir "${params.outdir}/08_functional_annotation/funannotate", mode: params.publish_mode

    input:
    path gff3
    path assembly
    path eggnog_anno
    path iprscan_xml
    path phobius_results

    output:
    path "final_anno", emit: anno_dir

    script:
    """
    funannotate annotate \\
      --gff ${gff3} \\
      --fasta ${assembly} \\
      -s ${params.species_name} \\
      --busco_db ${params.busco_lineage_meta} \\
      --eggnog ${eggnog_anno} \\
      --iprscan ${iprscan_xml} \\
      --phobius ${phobius_results} \\
      --cpus ${task.cpus} \\
      -o final_anno
    """

    stub:
    """
    mkdir -p final_anno
    echo "FUNANNOTATE_ANNOTATE stub completed" > final_anno/annotate.log
    """
}

// ============================================================
// Process: FUNANNOTATE_STATS
// Generate final annotation statistics
// ============================================================
process FUNANNOTATE_STATS {
    label 'functional_env'
    tag   "funannotate_stats"
    cpus  { 4 * task.attempt }
    memory { 8.GB * task.attempt }
    time   { 2.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::funannotate'
    container 'quay.io/biocontainers/funannotate:1.8.17--pyhdfd78af_5'
    publishDir "${params.outdir}/08_functional_annotation/funannotate", mode: params.publish_mode

    input:
    path assembly
    path gff3

    output:
    path "annotation_stats.json", emit: stats_json

    script:
    """
    funannotate util stats \\
      -f ${assembly} \\
      -g ${gff3} \\
      -o annotation_stats.json
    """

    stub:
    """
    touch annotation_stats.json
    echo "FUNANNOTATE_STATS stub completed"
    """
}

// ============================================================
// Process: FINAL_BUSCO
// Final BUSCO on annotated genome
// ============================================================
process FINAL_BUSCO {
    label 'functional_env'
    tag   "busco_annotated"
    cpus  { params.threads_busco * task.attempt }
    memory { 64.GB * task.attempt }
    time   { 48.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::busco'
    container 'quay.io/biocontainers/busco:6.1.0--pyhdfd78af_1'
    publishDir "${params.outdir}/08_functional_annotation/qc/busco", mode: params.publish_mode

    input:
    path gff3
    path assembly
    val  lineage

    output:
    path "busco_annotated_${lineage}/full_table.tsv", emit: busco_table

    script:
    """
    busco -i ${gff3} -m transcriptome -l ${lineage} -c ${task.cpus} \\
      --offline ${params.busco_offline_dir} \\
      -o busco_annotated_${lineage}
    """

    stub:
    """
    mkdir -p busco_annotated_${lineage}
    touch busco_annotated_${lineage}/full_table.tsv
    echo "FINAL_BUSCO stub completed"
    """
}

// // ============================================================
// // Workflow: functional_annotation
// // Phobius → InterProScan → eggNOG → funannotate annotate → stats
// // ============================================================
// workflow functional_annotation {

//     take:
//         assembly
//         gff3
//         proteins
//         funannotate_db

//     main:
//         PHOBIUS_RUN(proteins)
//         IPRSCAN_RUN(proteins)
//         EGGNOG_MAPPER(proteins, funannotate_db)
//         FUNANNOTATE_ANNOTATE(
//             gff3,
//             assembly,
//             EGGNOG_MAPPER.out.eggnog_anno,
//             IPRSCAN_RUN.out.iprscan_xml,
//             PHOBIUS_RUN.out.phobius_results
//         )
//         FUNANNOTATE_STATS(assembly, gff3)

//     emit:
//         anno_dir   = FUNANNOTATE_ANNOTATE.out.anno_dir
//         stats_json = FUNANNOTATE_STATS.out.stats_json
// }
