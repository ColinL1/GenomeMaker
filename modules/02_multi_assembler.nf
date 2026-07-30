// ============================================================
// 02_multi_assembler.nf
// Multi-assembler de novo assembly: NECAT / Flye / hifiasm / Raven
// ============================================================

nextflow.enable.dsl = 2

process ASSEMBLE_NECAT {
    label 'necat_env'
    tag   "necat"
    cpus  { params.resource_necat_cpus * task.attempt }
    memory { params.resource_necat_memory * task.attempt }
    time   { params.resource_necat_time * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    publishDir "${params.outdir}/02_assembly/necat", mode: params.publish_mode

    input:
    path assembly_reads
    path necat_config
    path necat_read_list

    output:
    path "assembly.fasta", emit: assembly

    script:
    """
    cp ${necat_config} necat_config.txt
    cp -r ${params.necat_path_bin} bin

    bin/necat.pl correct  necat_config.txt
    bin/necat.pl assemble necat_config.txt
    bin/necat.pl bridge   necat_config.txt

    cp \$(grep 'PROJECT=' necat_config.txt | sed 's/PROJECT=//g')/6-bridge_contigs/polished_contigs.fasta assembly.fasta
    """
    // script:
    // """
    // cp ${necat_config} necat_config.txt
    // # cp ${necat_read_list} necat_read_list.txt
    // export PATH=\$PATH:${params.necat_path_bin}

    // necat.pl correct  necat_config.txt
    // necat.pl assemble necat_config.txt
    // necat.pl bridge   necat_config.txt

    // # NECAT writes results under <PROJECT_NAME>/6-bridge_contigs/polished_contigs.fasta.
    // # PROJECT_NAME must come from params (must match the config's own PROJECT_NAME field).
    
    // cp \$(grep 'PROJECT=' necat_config.txt | sed 's/PROJECT=//g')/6-bridge_contigs/polished_contigs.fasta assembly.fasta
    // """

    stub:
    """
    touch assembly.fasta
    echo "ASSEMBLE_NECAT stub completed"
    """
}

process ASSEMBLE_FLYE {
    label 'necat_env'
    tag   "flye"
    cpus  { params.resource_flye_cpus * task.attempt }
    memory { params.resource_flye_memory * task.attempt }
    time   { params.resource_flye_time * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::flye'
    container 'quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1'
    publishDir "${params.outdir}/02_assembly/flye", mode: params.publish_mode

    input:
    path assembly_reads
    val  genome_size

    output:
    path "assembly.fasta", emit: assembly

    script:
    def genome_size_arg = genome_size ? "--genome-size ${genome_size}" : ""
    def iteration_args = params.iteration_flye ? "${params.iteration_flye}" : "10"
    """
    flye --nano-hq ${assembly_reads} \\
        ${genome_size_arg} \\
        --threads ${task.cpus} \\
        --out-dir flye_out \\
        --iterations ${iteration_args} \\
        --meta

    cp flye_out/assembly.fasta assembly.fasta
    """

    stub:
    """
    touch assembly.fasta
    echo "ASSEMBLE_FLYE stub completed"
    """
}

process ASSEMBLE_HIFIASM {
    label 'necat_env'
    tag   "hifiasm"
    cpus  { params.resource_hifiasm_cpus * task.attempt }
    memory { params.resource_hifiasm_memory * task.attempt }
    time   { params.resource_hifiasm_time * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::hifiasm'
    container 'quay.io/biocontainers/hifiasm:0.25.0--h5ca1c30_0'
    publishDir "${params.outdir}/02_assembly/hifiasm", mode: params.publish_mode

    input:
    path assembly_reads
    val  genome_size

    output:
    path "assembly.fasta", emit: assembly

    script:
    def genome_size_arg = genome_size ? "--hg-size ${(genome_size as double) / 1e6}m" : ""
    // -l2 (default purge level) rather than -l0: coral genomes are heterozygous
    // enough that disabling haplotig purging tends to retain both haplotypes
    // as separate contigs, inflating apparent size / BUSCO duplication.
    """
    mkdir -p hifiasm_out
    
    hifiasm -o hifiasm_out/asm -t ${task.cpus} --ont -l2 ${genome_size_arg} ${assembly_reads}

    awk '/^S/{print ">"\$2"\\n"\$3}' hifiasm_out/asm.bp.p_ctg.gfa > assembly.fasta
    """

    stub:
    """
    touch assembly.fasta
    echo "ASSEMBLE_HIFIASM stub completed"
    """
}

process ASSEMBLE_RAVEN {
    label 'necat_env'
    tag   "raven"
    cpus  { params.resource_raven_cpus * task.attempt }
    memory { params.resource_raven_memory * task.attempt }
    time   { params.resource_raven_time * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::raven'
    container 'quay.io/biocontainers/raven-assembler:1.8.3--h5ca1c30_3'
    publishDir "${params.outdir}/02_assembly/raven", mode: params.publish_mode

    input:
    path assembly_reads

    output:
    path "assembly.fasta", emit: assembly

    script:
    """
    raven -t ${task.cpus} ${assembly_reads} > assembly.fasta
    """

    stub:
    """
    touch assembly.fasta
    echo "ASSEMBLE_RAVEN stub completed"
    """
}

// ============================================================
// QC processes. Filenames embed ${name} so multiple assemblers'
// reports can be safely staged together into COMBINE_COMPARISON
// and RANK_ASSEMBLIES without collisions.
// ============================================================

process BUSCO_ASSEMBLY {
    label 'necat_env'
    tag   "busco_${name}_${lineage}"
    cpus  { params.threads_busco * task.attempt }
    memory { 64.GB * task.attempt }
    time   { 48.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::busco'
    container 'quay.io/biocontainers/busco:6.1.0--pyhdfd78af_1'
    publishDir "${params.outdir}/02_assembly/qc/busco", mode: params.publish_mode

    input:
    tuple val(name), path(assembly)
    val  lineage

    output:
    tuple val(name), path("${name}.${lineage}.full_table.tsv"),    emit: busco_table
    tuple val(name), path("${name}.${lineage}.short_summary.txt"), emit: busco_summary
    tuple val(name), path("${name}.${lineage}.short_summary.json"), emit: busco_summary_json

    script:
    def busco_dir = params.busco_offline_dir ? "--download_path  ${params.busco_offline_dir}" : ""
    """
    busco -i ${assembly} -m geno -l ${lineage} -c ${task.cpus} -o busco_output ${busco_dir}
    cp busco_output/run_${lineage}/full_table.tsv ${name}.${lineage}.full_table.tsv
    cp busco_output/short_summary.specific.${lineage}.busco_output.txt ${name}.${lineage}.short_summary.txt
    cp busco_output/short_summary.specific.${lineage}.busco_output.json ${name}.${lineage}.short_summary.json
    """

    stub:
    """
    touch ${name}.${lineage}.full_table.tsv ${name}.${lineage}.short_summary.txt ${name}.${lineage}.short_summary.json
    echo "BUSCO_ASSEMBLY stub completed"
    """
}

process GFASTATS_QC {
    label 'necat_env'
    tag   "gfastats_${name}"
    cpus  { 4 * task.attempt }
    memory { 8.GB * task.attempt }
    time   { 2.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::gfastats'
    container 'quay.io/biocontainers/gfastats:1.3.11--h077b44d_0'
    publishDir "${params.outdir}/02_assembly/qc/gfastats", mode: params.publish_mode

    input:
    tuple val(name), path(assembly)

    output:
    tuple val(name), path("${name}.gfastats.txt"), emit: gfastats_report

    script:
    """
    gfastats ${assembly} > ${name}.gfastats.txt
    """

    stub:
    """
    touch ${name}.gfastats.txt
    echo "GFASTATS_QC stub completed"
    """
}

process QUAST_QC {
    label 'necat_env'
    tag   "quast_${name}"
    cpus  { 4 * task.attempt }
    memory { 16.GB * task.attempt }
    time   { 4.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::quast'
    container 'quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2'
    publishDir "${params.outdir}/02_assembly/qc/quast", mode: params.publish_mode

    input:
    tuple val(name), path(assembly)

    output:
    tuple val(name), path("${name}_quast_report"), emit: quast_report

    script:
    """
    quast.py ${assembly} -o ${name}_quast_report --threads ${task.cpus} --min-contig 1000
    """

    stub:
    """
    mkdir -p ${name}_quast_report
    touch ${name}_quast_report/report.tsv
    echo "QUAST_QC stub completed"
    """
}

// ============================================================
// Process: COMBINE_COMPARISON
// Aggregate per-assembler QC into comparison reports.
// Parses BUSCO short_summary.txt (the "C:91.4%[S:...,D:...],F:...,M:...,n:..."
// line), NOT full_table.tsv -- full_table.tsv has no such summary line,
// so the previous version's parser silently matched nothing.
// ============================================================
process COMBINE_COMPARISON {
    tag   "combine_comparison"
    cpus  { 2 * task.attempt }
    memory { 4.GB * task.attempt }
    time   { 1.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'conda-forge::pandas conda-forge::jinja2'
    container 'quay.io/biocontainers/mulled-v2-1c6be8ad49e4dfe8ab70558e8fb200d7b2fd7509:5900b4e68c4051137fffd99165b00e98f810acae-0'
    publishDir "${params.outdir}/02_assembly/comparison", mode: params.publish_mode

    input:
    path gfastats_reports    // ${name}.gfastats.txt files
    path busco_summaries     // ${name}.${lineage}.short_summary.txt files (host lineage)

    output:
    path "assembly_comparison.tsv",  emit: comparison_tsv
    path "assembly_comparison.json", emit: comparison_json

    script:
    """
    python3 ${params.bin_path}combine_comparison.py
    """

    stub:
    """
    touch assembly_comparison.tsv assembly_comparison.json
    echo "COMBINE_COMPARISON stub completed"
    """
}
