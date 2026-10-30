// ============================================================
// 00_preprocess_reads.nf
// ONT read preprocessing: porechop adapter removal, QC, chopper split
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: PORECHOP
// Remove adapters from ONT reads
// ============================================================
process PORECHOP {
    label 'reads_env'
    tag   "porechop_${reads.simpleName}"

    conda 'bioconda::porechop'
    container 'quay.io/biocontainers/porechop:0.2.4--py311he264feb_9'
    publishDir "${params.outdir}/00_preprocessing/porechop", mode: params.publish_mode_intermediate

    input:
    path reads

    output:
    path "ONT_reads_pc.fastq.gz", emit: reads_pc
    path "versions.yml", emit: versions

    script:
    """
    porechop \\
      --input ${reads} \\
      -o ONT_reads_pc.fastq.gz \\
      --discard_middle
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        porechop: "\$(porechop --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """
    stub:
    """
    touch ONT_reads_pc.fastq.gz
    echo "PORECHOP stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        porechop: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: READ_QC
// Quality control: FastQC, MultiQC, NanoPlot
// ============================================================
process FASTQC {
    label 'reads_env'
    label 'res_read_tool'
    tag   { "read_qc_${stage}" }

    conda 'bioconda::fastqc'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir { "${params.outdir}/00_preprocessing/qc/${stage}/fastqc" }, mode: params.publish_mode

    input:
    path reads_pc
    val  stage

    output:
    path "qc_fastqc_${stage}", emit: qc_dir
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p qc_fastqc_${stage}
    fastqc -t ${task.cpus} ${reads_pc} -o qc_fastqc_${stage}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: "\$(fastqc --version 2>&1 | sed 's/FastQC v//' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """
    stub:
    """
    mkdir -p qc_fastqc_${stage}
    touch qc_fastqc_${stage}/${reads_pc}_fastqc.zip
    touch qc_fastqc_${stage}/${reads_pc}_fastqc.html
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: "stub"
    END_VERSIONS
    """
}

process NANOPLOT {
    label 'reads_env'
    label 'res_read_tool'
    tag   { "read_qc_${stage}" }

    conda 'nanoplot'
    container 'quay.io/biocontainers/nanoplot:1.47.1--pyhdfd78af_0'
    publishDir { "${params.outdir}/00_preprocessing/qc/${stage}/nanoplot" }, mode: params.publish_mode

    input:
    path reads_pc
    val  stage

    output:
    path "qc_nanoplot_${stage}", emit: qc_dir
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p qc_nanoplot_${stage}
    NanoPlot --fastq ${reads_pc} --plots dot --legacy hex --N50 -o qc_nanoplot_${stage} -t ${task.cpus}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: "\$(NanoPlot --version 2>&1 | sed 's/NanoPlot //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """
    stub:
    """
    mkdir -p qc_nanoplot_${stage}
    touch qc_nanoplot_${stage}/NanoPlot-report.html
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: "stub"
    END_VERSIONS
    """
}

process NANOSTAT {
    label 'reads_env'
    label 'res_read_tool'
    tag   { "read_qc_${stage}" }

    conda 'nanostat'
    container 'quay.io/biocontainers/nanostat:0.1.5--py35_0'
    publishDir { "${params.outdir}/00_preprocessing/qc/${stage}/nanostat" }, mode: params.publish_mode

    input:
    path reads_pc
    val  stage

    output:
    path "qc_nanostat_${stage}", emit: qc_dir
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p qc_nanostat_${stage}
    NanoStat --fastq ${reads_pc} -t ${task.cpus} > qc_nanostat_${stage}/nanostat_${stage}.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanostat: "\$(NanoStat --version 2>&1 | sed 's/NanoStat //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """
    stub:
    """
    mkdir -p qc_nanostat_${stage}
    touch qc_nanostat_${stage}/nanostat_${stage}.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanostat: "stub"
    END_VERSIONS
    """
}

process MULTIQC {
    label 'reads_env'
    label 'res_read_tool'
    tag   { "read_qc_${stage}" }

    conda 'bioconda::multiqc'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir { "${params.outdir}/00_preprocessing/qc/${stage}/multiqc" }, mode: params.publish_mode

    input:
    path fastqc_dirs,   stageAs: 'fastqc_input_*'
    path nanoplot_dirs, stageAs: 'nanoplot_input_*'
    path nanostat_dirs, stageAs: 'nanostat_input_*'
    val  stage

    output:
    path "qc", emit: qc_dir
    path "versions.yml", emit: versions

    script:
    """
    multiqc fastqc_input_* nanoplot_input_* nanostat_input_* -o qc --interactive --force
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: "\$(multiqc --version 2>&1 | sed 's/.*version //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """
    stub:
    """
    mkdir -p qc
    touch qc/multiqc_report.${stage}.testrun.html
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: CHOPPER_HQ_READS
// High-quality read subset used ONLY for k-mer genome size estimation.
//
// WHY A SEPARATE PROCESS, not a third CHOPPER_SPLIT output: CHOPPER_SPLIT's
// outputs feed purge_dups, the mitogenome stage and all of contamination
// filtering, so touching its script invalidates the cache for most of the
// pipeline. Deriving the k-mer read set here keeps that cache intact.
//
// WHY FILTER AT ALL: GenomeScope models a k-mer spectrum assuming most k-mers are
// either genomic or clearly erroneous. Raw ONT error rates smear the error peak
// into the homozygous peak, which inflates the apparent genome size. Restricting
// to higher-quality reads sharpens both peaks. The estimate remains advisory for
// a long-read-only dataset -- say so in methods.
// ============================================================
process CHOPPER_HQ_READS {
    label 'reads_env'
    tag   "chopper_hq_q${params.kmer_min_quality}"

    conda 'bioconda::chopper'
    container 'quay.io/biocontainers/chopper:0.13.0--h7f49ad2_0'
    publishDir "${params.outdir}/00_preprocessing/reads", mode: params.publish_mode_intermediate

    input:
    path reads_pc

    output:
    path "kmer_reads_q${params.kmer_min_quality}.fastq.gz", emit: hq_reads
    path "versions.yml", emit: versions

    script:
    def half = Math.max(1, (task.cpus / 2) as int)
    """
    (
      set -o pipefail
      gunzip -c ${reads_pc} \\
        | chopper -q ${params.kmer_min_quality} --threads ${half} \\
        | gzip > kmer_reads_q${params.kmer_min_quality}.fastq.gz
    )

    if [ ! -s kmer_reads_q${params.kmer_min_quality}.fastq.gz ]; then
        echo "ERROR: no reads survived the q${params.kmer_min_quality} filter -- lower params.kmer_min_quality" >&2
        exit 1
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: "\$(chopper --version 2>&1 | sed 's/chopper //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch kmer_reads_q${params.kmer_min_quality}.fastq.gz
    echo "CHOPPER_HQ_READS stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: CHOPPER_SPLIT
// Split reads into an assembly set (q3, >=1 kb) and a polishing set (q5, >=500 bp)
// Optionally filter symbiont contamination
// ============================================================
process CHOPPER_SPLIT {
    label 'reads_env'
    tag   "chopper_split"

    conda 'bioconda::chopper'
    container 'quay.io/biocontainers/chopper:0.13.0--h7f49ad2_0'
    publishDir "${params.outdir}/00_preprocessing/reads", mode: params.publish_mode_intermediate

    input:
    path reads_pc
    path symbiont_ref  // optional; pass [] when not provided

    output:
    path "assembly_reads.fastq.gz",  emit: assembly_reads
    path "polishing_reads.fastq.gz", emit: polishing_reads
    path "versions.yml", emit: versions

    script:
    def contam_opt = !params.run_symbiont_filter.toString().toBoolean() || !symbiont_ref
        ? ''
        : "--contam ${symbiont_ref}"
    def half_cpus = Math.max(1, (task.cpus / 2) as int)
    """
    # ONE decompression, tee'd into BOTH chopper streams, which run concurrently.
    # The assembly stream reads from a FIFO in the background, so its exit status
    # can be collected with `wait PID` (a >(...) process substitution cannot be).
    #
    # gzip -2 rather than the default -6: these are intermediate read sets that
    # exist only to feed the assembler and the polisher, and this container ships
    # no pigz/bgzip, so compression is single-threaded per stream.
    set -o pipefail
    mkfifo assembly_stream
    chopper -q 3 -l 1000 ${contam_opt} --threads ${half_cpus} < assembly_stream \\
        | gzip -2 > assembly_reads.fastq.gz &
    assembly_pid=\$!

    gunzip -c ${reads_pc} \\
        | tee assembly_stream \\
        | chopper -q 5 -l 500 ${contam_opt} --threads ${half_cpus} \\
        | gzip -2 > polishing_reads.fastq.gz

    wait \$assembly_pid
    rm -f assembly_stream

    for f in assembly_reads.fastq.gz polishing_reads.fastq.gz; do
        if [ ! -s "\$f" ]; then
            echo "ERROR: \$f is empty after chopper -- check the quality/length thresholds" >&2
            exit 1
        fi
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: "\$(chopper --version 2>&1 | sed 's/chopper //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch assembly_reads.fastq.gz polishing_reads.fastq.gz
    echo "CHOPPER_SPLIT stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: "stub"
    END_VERSIONS
    """
}
