// ============================================================
// 04_assembly_polish.nf
// Assembly polishing: Racon (ONT) → Medaka (ONT) → funannotate clean/sort
// Final BUSCO on polished assembly
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: RACON_POLISH
// Iterative Racon polishing with ONT reads
// ============================================================
process RACON_POLISH {
    label 'polish_env'
    tag   "racon_polish"
    cpus  { params.threads_default * task.attempt }
    memory { 64.GB * task.attempt }
    time   { 12.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::racon minimap2'
    container 'quay.io/biocontainers/racon:1.5.0--h077b44d_8'
    publishDir "${params.outdir}/04_polishing/racon", mode: params.publish_mode

    input:
    path assembly
    path reads

    output:
    path "polished.fasta", emit: polished_assembly

    script:
    """
    # Map reads to assembly
    minimap2 -a -x map-ont -t ${task.cpus} \\
      ${assembly} ${reads} \\
      -o aligned.sam

    # Racon polishing (1 iteration, can be increased)
    racon -u -m 3 -x -5 -g -4 -w 500 -t ${task.cpus} \\
      ${reads} aligned.sam ${assembly} \\
      > polished.fasta
    """

    stub:
    """
    touch polished.fasta
    echo "RACON_POLISH stub completed"
    """
}

// ============================================================
// Process: MEDAKA_POLISH
// Medaka consensus polishing with ONT reads
// ============================================================
process MEDAKA_POLISH {
    label 'polish_env'
    tag   "medaka_polish"
    cpus  { params.threads_default * task.attempt }
    memory { 32.GB * task.attempt }
    time   { 12.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::medaka minimap2'
    container 'quay.io/biocontainers/medaka:2.2.2--py312h3050eb1_0'
    publishDir "${params.outdir}/04_polishing/medaka", mode: params.publish_mode

    input:
    path assembly
    path reads
    val  medaka_model   // e.g. 'r1041_e82_400bps_hac_v4.3.0', auto-selected if null

    output:
    path "medaka_polished.fasta", emit: medaka_assembly

    script:
    def model = medaka_model ?: 'auto'

    """
    # Map reads to assembly
    minimap2 -a -x map-ont -t ${task.cpus} \\
      ${assembly} ${reads} \\
      -o aligned.sam

    # Convert to sorted BAM
    samtools view -Sb -@ ${task.cpus} aligned.sam | \\
      samtools sort -@ ${task.cpus} -o aligned.sorted.bam
    samtools index aligned.sorted.bam

    # Medaka consensus
    if [ "${model}" = "auto" ]; then
      # Auto-detect model based on basecaller (simplified)
      model=\$(medaka list_models | head -1 | awk '{print \$1}')
      echo "Auto-selected medaka model: \${model}"
    fi

    medaka.consensus -i aligned.sorted.bam -d ${assembly} \\
      -m ${model} -t ${task.cpus} -o medaka_out

    cp medaka_out/consensus.fasta medaka_polished.fasta
    """

    stub:
    """
    touch medaka_polished.fasta
    echo "MEDAKA_POLISH stub completed"
    """
}

// ============================================================
// Process: FUNANNOTATE_CLEAN_SORT
// Clean (remove short contigs) and sort assembly
// ============================================================
process FUNANNOTATE_CLEAN_SORT {
    label 'polish_env'
    tag   "funannotate_clean_sort"
    cpus  { 4 * task.attempt }
    memory { 16.GB * task.attempt }
    time   { 4.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::funannotate'
    container 'quay.io/biocontainers/funannotate:1.8.17--pyhdfd78af_5'
    publishDir "${params.outdir}/04_polishing/funannotate", mode: params.publish_mode

    input:
    path assembly
    val  species_name

    output:
    path "clean_sorted.fasta", emit: final_assembly

    script:
    """
    # Clean: remove contigs shorter than 200 bp
    funannotate clean -i ${assembly} -m 200 -o clean.fasta

    # Sort: rename and sort scaffolds by size (largest first)
    funannotate sort -i clean.fasta -o clean_sorted.fasta -b ${species_name}_scaffold
    """

    stub:
    """
    touch clean_sorted.fasta
    echo "FUNANNOTATE_CLEAN_SORT stub completed"
    """
}

// ============================================================
// Process: BUSCO_FINAL
// Final BUSCO on polished assembly
// ============================================================
process BUSCO_FINAL {
    label 'polish_env'
    tag   "busco_final_${lineage}"
    cpus  { params.threads_busco * task.attempt }
    memory { 64.GB * task.attempt }
    time   { 48.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::busco'
    container 'quay.io/biocontainers/busco:6.1.0--pyhdfd78af_1'
    publishDir "${params.outdir}/04_polishing/qc/busco", mode: params.publish_mode

    input:
    path assembly
    val  lineage

    output:
    path "busco_final_${lineage}/full_table.tsv", emit: busco_table
    path "busco_final_${lineage}/short_summary.txt", emit: busco_summary

    script:
    """
    busco -i ${assembly} -m geno -l ${lineage} -c ${task.cpus} \\
      -o busco_final_${lineage}
    """

    stub:
    """
    mkdir -p busco_final_${lineage}
    touch busco_final_${lineage}/full_table.tsv busco_final_${lineage}/short_summary.txt
    echo "BUSCO_FINAL stub completed"
    """
}

process BUSCO_FINAL_META {
    label 'polish_env'
    tag   "busco_final_meta_${lineage}"
    cpus  { params.threads_busco * task.attempt }
    memory { 64.GB * task.attempt }
    time   { 48.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::busco'
    container 'quay.io/biocontainers/busco:6.1.0--pyhdfd78af_1'
    publishDir "${params.outdir}/04_polishing/qc/busco", mode: params.publish_mode

    input:
    path assembly
    val  lineage

    output:
    path "busco_final_meta_${lineage}/full_table.tsv", emit: busco_table
    path "busco_final_meta_${lineage}/short_summary.txt", emit: busco_summary

    script:
    """
    busco -i ${assembly} -m geno -l ${lineage} -c ${task.cpus} \\
      -o busco_final_meta_${lineage}
    """

    stub:
    """
    mkdir -p busco_final_meta_${lineage}
    touch busco_final_meta_${lineage}/full_table.tsv busco_final_meta_${lineage}/short_summary.txt
    echo "BUSCO_FINAL_META stub completed"
    """
}

// // ============================================================
// // Workflow: assembly_polishing
// // Racon → Medaka → clean/sort → final BUSCO
// // ============================================================
// workflow assembly_polishing {

//     take:
//         assembly
//         reads
//         busco_lineage_euk
//         busco_lineage_meta

//     main:
//         RACON_POLISH(assembly, reads)
//         MEDAKA_POLISH(RACON_POLISH.out.polished_assembly, reads, channel.value([]))
//         FUNANNOTATE_CLEAN_SORT(MEDAKA_POLISH.out.medaka_assembly, params.species_name)
//         BUSCO_FINAL(FUNANNOTATE_CLEAN_SORT.out.final_assembly, busco_lineage_euk)
//         BUSCO_FINAL_META(FUNANNOTATE_CLEAN_SORT.out.final_assembly, busco_lineage_meta)

//     emit:
//         final_assembly  = FUNANNOTATE_CLEAN_SORT.out.final_assembly
//         busco_tables    = BUSCO_FINAL.out.busco_table.mix(BUSCO_FINAL_META.out.busco_table)
//         busco_summaries = BUSCO_FINAL.out.busco_summary.mix(BUSCO_FINAL_META.out.busco_summary)
// }
