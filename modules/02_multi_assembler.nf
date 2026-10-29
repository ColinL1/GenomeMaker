// ============================================================
// 02_multi_assembler.nf
// Multi-assembler de novo assembly: NECAT / Flye / hifiasm / Raven
// ============================================================

nextflow.enable.dsl = 2

process ASSEMBLE_NECAT {
    label 'assembly_env'
    tag   "necat"

    publishDir "${params.outdir}/02_assembly/necat", mode: params.publish_mode

    input:
    path assembly_reads
    path necat_config
    path necat_read_list
    path necat_bin_dir, stageAs: 'necat_bin_src'

    output:
    path "assembly.fasta", emit: assembly
    path "versions.yml", emit: versions

    script:
    """
    cp ${necat_config} necat_config.txt
    cp -r ${necat_bin_dir} bin

    bin/necat.pl correct  necat_config.txt
    bin/necat.pl assemble necat_config.txt
    bin/necat.pl bridge   necat_config.txt

    cp \$(grep 'PROJECT=' necat_config.txt | sed 's/PROJECT=//g')/6-bridge_contigs/polished_contigs.fasta assembly.fasta
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        necat: "\$(bin/necat.pl 2>&1 | grep -im1 version | sed 's/.*ersion *//' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch assembly.fasta
    echo "ASSEMBLE_NECAT stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        necat: "stub"
    END_VERSIONS
    """
}

process ASSEMBLE_FLYE {
    label 'assembly_env'
    tag   "flye"

    conda 'bioconda::flye'
    container 'quay.io/biocontainers/flye:2.9.6--py313h7fbb527_1'
    publishDir "${params.outdir}/02_assembly/flye", mode: params.publish_mode

    input:
    path assembly_reads
    val  genome_size

    output:
    path "assembly.fasta", emit: assembly
    path "versions.yml", emit: versions

    script:
    def genome_size_arg = genome_size ? "--genome-size ${genome_size}" : ""
    // --meta is metagenome mode: it relaxes Flye's uniform-coverage assumption,
    // which suits a holobiont sample carrying Symbiodiniaceae and bacteria at
    // very different depths. It also changes repeat resolution and interacts with
    // --genome-size, so it is a deliberate, documented choice rather than a
    // default -- state it in methods. Toggle it off for a clean host-only library.
    def meta_arg = params.flye_meta.toString().toBoolean() ? "--meta" : ""
    """
    flye --nano-hq ${assembly_reads} \\
        ${genome_size_arg} \\
        --threads ${task.cpus} \\
        --out-dir flye_out \\
        --iterations ${params.iteration_flye} \\
        ${meta_arg}

    cp flye_out/assembly.fasta assembly.fasta
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: "\$(flye --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch assembly.fasta
    echo "ASSEMBLE_FLYE stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: "stub"
    END_VERSIONS
    """
}

process ASSEMBLE_HIFIASM {
    label 'assembly_env'
    tag   "hifiasm"

    conda 'bioconda::hifiasm'
    container 'quay.io/biocontainers/hifiasm:0.25.0--h5ca1c30_0'
    publishDir "${params.outdir}/02_assembly/hifiasm", mode: params.publish_mode

    input:
    path assembly_reads
    val  genome_size

    output:
    path "assembly.fasta", emit: assembly
    path "versions.yml", emit: versions

    script:
    def genome_size_arg = genome_size ? "--hg-size ${(genome_size as double) / 1e6}m" : ""
    // -l2 (default purge level) rather than -l0: coral genomes are heterozygous
    // enough that disabling haplotig purging tends to retain both haplotypes
    // as separate contigs, inflating apparent size / BUSCO duplication.
    """
    mkdir -p hifiasm_out
    
    hifiasm -o hifiasm_out/asm -t ${task.cpus} --ont -l2 ${genome_size_arg} ${assembly_reads}

    awk '/^S/{print ">"\$2"\\n"\$3}' hifiasm_out/asm.bp.p_ctg.gfa > assembly.fasta
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hifiasm: "\$(hifiasm --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch assembly.fasta
    echo "ASSEMBLE_HIFIASM stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hifiasm: "stub"
    END_VERSIONS
    """
}

process ASSEMBLE_RAVEN {
    label 'assembly_env'
    tag   "raven"

    conda 'bioconda::raven'
    container 'quay.io/biocontainers/raven-assembler:1.8.3--h5ca1c30_3'
    publishDir "${params.outdir}/02_assembly/raven", mode: params.publish_mode

    input:
    path assembly_reads

    output:
    path "assembly.fasta", emit: assembly
    path "versions.yml", emit: versions

    script:
    """
    raven -t ${task.cpus} ${assembly_reads} > assembly.fasta
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        raven: "\$(raven --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch assembly.fasta
    echo "ASSEMBLE_RAVEN stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        raven: "stub"
    END_VERSIONS
    """
}

// ============================================================
// QC processes. Filenames embed ${name} so multiple assemblers'
// reports can be safely staged together into COMBINE_COMPARISON
// and RANK_ASSEMBLIES without collisions.
// ============================================================

process BUSCO_ASSEMBLY {
    label 'assembly_env'
    tag   "busco_${name}_${lineage}"

    conda 'bioconda::busco'
    container 'quay.io/biocontainers/busco:6.1.0--pyhdfd78af_1'
    publishDir "${params.outdir}/02_assembly/qc/busco", mode: params.publish_mode

    input:
    tuple val(name), path(assembly)
    val  lineage
    path busco_db_dir

    output:
    tuple val(name), path("${name}.${lineage}.full_table.tsv"),    emit: busco_table
    tuple val(name), path("${name}.${lineage}.short_summary.txt"), emit: busco_summary
    tuple val(name), path("${name}.${lineage}.short_summary.json"), emit: busco_summary_json
    path "versions.yml", emit: versions

    script:
    // --offline is required alongside --download_path, not optional. Without it
    // BUSCO still contacts busco-data.ezlab.org to version-check the lineage and
    // errors out (it does not fall back gracefully) on a node with no egress.
    // It also pins the dataset version: without --offline a silent upstream update
    // mid-project makes BUSCO scores non-comparable across runs of this pipeline.
    def busco_dir = params.busco_offline_dir ? "--offline --download_path ${busco_db_dir}" : ""
    def busco_mem = "${task.memory.toGiga()}g"
    """
    export _JAVA_OPTIONS="-Xmx${busco_mem}"
    busco -i ${assembly} -m geno -l ${lineage} -c ${task.cpus} -o busco_output ${busco_dir}
    cp busco_output/run_${lineage}/full_table.tsv ${name}.${lineage}.full_table.tsv
    cp busco_output/short_summary.specific.${lineage}.busco_output.txt ${name}.${lineage}.short_summary.txt
    cp busco_output/short_summary.specific.${lineage}.busco_output.json ${name}.${lineage}.short_summary.json
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: "\$(busco --version 2>&1 | sed 's/BUSCO //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch ${name}.${lineage}.full_table.tsv ${name}.${lineage}.short_summary.txt ${name}.${lineage}.short_summary.json
    echo "BUSCO_ASSEMBLY stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: "stub"
    END_VERSIONS
    """
}

process GFASTATS_QC {
    label 'assembly_env'
    label 'res_light'
    tag   "gfastats_${name}"

    conda 'bioconda::gfastats'
    container 'quay.io/biocontainers/gfastats:1.3.11--h077b44d_0'
    publishDir "${params.outdir}/02_assembly/qc/gfastats", mode: params.publish_mode

    input:
    tuple val(name), path(assembly)

    output:
    tuple val(name), path("${name}.gfastats.txt"), emit: gfastats_report
    path "versions.yml", emit: versions

    script:
    """
    gfastats ${assembly} > ${name}.gfastats.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gfastats: "\$(gfastats --version 2>&1 | sed 's/.*v//' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch ${name}.gfastats.txt
    echo "GFASTATS_QC stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gfastats: "stub"
    END_VERSIONS
    """
}

process QUAST_QC {
    label 'assembly_env'
    label 'res_mid_4h'
    tag   "quast_${name}"

    conda 'bioconda::quast'
    container 'quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2'
    publishDir "${params.outdir}/02_assembly/qc/quast", mode: params.publish_mode

    input:
    tuple val(name), path(assembly)

    output:
    tuple val(name), path("${name}_quast_report"), emit: quast_report
    path "versions.yml", emit: versions

    script:
    """
    quast.py ${assembly} -o ${name}_quast_report --threads ${task.cpus} --min-contig 1000
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: "\$(quast.py --version 2>&1 | grep -m1 -o 'QUAST v[0-9.]*' | sed 's/QUAST v//' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${name}_quast_report
    touch ${name}_quast_report/report.tsv
    echo "QUAST_QC stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: "stub"
    END_VERSIONS
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
    label 'res_tiny'
    tag   "combine_comparison"

    conda 'conda-forge::pandas conda-forge::jinja2'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/02_assembly/comparison", mode: params.publish_mode

    input:
    path gfastats_reports    // ${name}.gfastats.txt files
    path busco_summaries     // ${name}.${lineage}.short_summary.txt files (host lineage)

    output:
    path "assembly_comparison.tsv",  emit: comparison_tsv
    path "assembly_comparison.json", emit: comparison_json
    path "versions.yml", emit: versions

    script:
    """
    combine_comparison.py
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "\$(python3 --version 2>&1 | sed 's/Python //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch assembly_comparison.tsv assembly_comparison.json
    echo "COMBINE_COMPARISON stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
