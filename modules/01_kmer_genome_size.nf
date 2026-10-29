// ============================================================
// 01_kmer_genome_size.nf
// K-mer counting, histogram, and GenomeScope 2.0 estimation
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: MERYL_COUNT
// Count k-mers from reads
// ============================================================
process MERYL_COUNT {
    label 'kmer_env'
    tag   "meryl_count"

    conda 'bioconda::meryl=1.4.2'
    container 'quay.io/biocontainers/meryl:1.4.2--hd19868c_0'
    publishDir "${params.outdir}/01_kmer_genome_size/meryl", mode: params.publish_mode_intermediate

    input:
    path reads_pc

    output:
    path "kmer_db.meryl", emit: meryl_db
    path "versions.yml", emit: versions

    script:
    // meryl defaults to ALL cores and ALL memory on the host, ignoring the
    // scheduler's allocation -- under SLURM that oversubscribes the node and can
    // get the job killed. threads=/memory= bind it to what Nextflow asked for.
    // memory= is in GB.
    """
    meryl count \\
      k=${params.kmer_size} \\
      threads=${task.cpus} \\
      memory=${task.memory.toGiga()} \\
      output kmer_db.meryl \\
      ${reads_pc}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: "\$(meryl --version 2>&1 | sed 's/meryl //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p kmer_db.meryl
    echo "MERYL_COUNT stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: MERYL_HIST
// Generate k-mer frequency histogram
// ============================================================
process MERYL_HIST {
    label 'kmer_env'
    label 'res_mid_6h'
    tag   "meryl_histogram"

    conda 'bioconda::meryl=1.4.2'
    container 'quay.io/biocontainers/meryl:1.4.2--hd19868c_0'
    publishDir "${params.outdir}/01_kmer_genome_size/meryl", mode: params.publish_mode_intermediate

    input:
    path meryl_db

    output:
    path "kmer_db.meryl.hist", emit: meryl_hist
    path "versions.yml", emit: versions

    script:
    """
    meryl histogram threads=${task.cpus} kmer_db.meryl > kmer_db.meryl.hist
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: "\$(meryl --version 2>&1 | sed 's/meryl //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch kmer_db.meryl.hist
    echo "MERYL_HIST stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: GENOMESCOPE
// Run GenomeScope 2.0 to estimate genome properties
// ============================================================
process GENOMESCOPE {
    label 'kmer_env'
    label 'res_small_2h'
    tag   "genomescope"

    conda 'bioconda::genomescope2=2.1.0'
    container 'quay.io/biocontainers/genomescope2:2.1.0--py313r44hdfd78af_0'
    publishDir "${params.outdir}/01_kmer_genome_size/genomescope", mode: params.publish_mode

    input:
    path meryl_hist

    output:
    path "genomescope_output_${params.kmer_size}", emit: genomescope_out
    path "versions.yml", emit: versions

    script:
    // GenomeScope2 from a container, replacing the vendored config/bin/genomescope.R
    // -- that script was GenomeScope *1.0* (its own header cites
    // schatzlab/genomescope) despite every comment here claiming 2.0. v1 has no
    // ploidy model, which matters for a highly heterozygous coral.
    //
    // The old positional call also passed a read length of 150 for ONT reads that
    // are kilobases long. GenomeScope2 takes no read-length argument at all -- it
    // only ever affected coverage reporting, never the haploid-length estimate --
    // so the parameter disappears rather than needing a corrected value.
    """
    awk 'BEGIN{OFS=" "}{print \$1, \$2}' ${meryl_hist} > kmer_db.meryl.hist.clean

    genomescope2 \\
      -i kmer_db.meryl.hist.clean \\
      -o genomescope_output_${params.kmer_size} \\
      -k ${params.kmer_size} \\
      -p ${params.ploidy}

    if [ ! -s genomescope_output_${params.kmer_size}/summary.txt ]; then
        echo "ERROR: GenomeScope2 produced no summary.txt -- the model did not converge." >&2
        echo "       Check progress.txt in that directory and the k-mer histogram." >&2
        exit 1
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomescope2: "\$(genomescope2 --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """
    stub:
    """
    mkdir -p genomescope_output_${params.kmer_size}
    printf 'Genome Haploid Length   600,000,000 bp   610,000,000 bp\\n' > genomescope_output_${params.kmer_size}/summary.txt
    echo "GENOMESCOPE stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomescope2: "stub"
    END_VERSIONS
    """
}

process PARSE_GENOME_SIZE {
    tag "parse_genomescope"
    label 'kmer_env'

    input:
    path genomescope_out

    output:
    env 'GENOME_SIZE', emit: genome_size
    path "versions.yml", emit: versions

    script:
    """
    SUMMARY=\$(ls ${genomescope_out}/*summary.txt | head -n1)

    GENOME_SIZE=\$(grep -m1 'Genome Haploid Length' "\$SUMMARY" \\
        | sed -E 's/^Genome Haploid Length[[:space:]]+([0-9,]+) bp.*/\\1/' \\
        | tr -d ',')

    if [ -z "\$GENOME_SIZE" ]; then
        echo "ERROR: could not parse Genome Haploid Length from \$SUMMARY" >&2
        exit 1
    fi
    echo "Parsed genome size: \$GENOME_SIZE bp"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: "\$(bash --version 2>&1 | head -1 | sed 's/.*version //;s/ .*//' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    GENOME_SIZE=600000000
    echo "PARSE_GENOME_SIZE stub completed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: "stub"
    END_VERSIONS
    """
}
