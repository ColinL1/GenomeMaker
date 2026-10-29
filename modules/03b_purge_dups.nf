// ============================================================
// 03b_purge_dups.nf
// Standard purge_dups pipeline (read-depth + self-alignment based
// haplotig purging), using ONT coverage. Runs on the HOST bin of the
// BlobToolKit split, so the coverage cutoffs come from host contigs only.
// Then the purge QC gate and the length filter (>= 200 bp).
// ============================================================

nextflow.enable.dsl = 2

process PURGE_DUPS_COVERAGE {
    label 'purge_dups_env'
    label 'res_purge_dups'
    tag   "purge_dups_coverage"

    conda 'bioconda::purge_dups=1.2.6 bioconda::minimap2'
    container 'quay.io/biocontainers/purge_dups:1.2.6--h577a1d6_3'
    publishDir "${params.outdir}/03b_purge_dups", mode: params.publish_mode_intermediate

    input:
    path assembly
    path reads   // ONT reads (assembly_reads)

    output:
    path "PB.base.cov", emit: base_cov
    path "PB.stat",     emit: stat
    path "cutoffs",     emit: cutoffs
    path "versions.yml", emit: versions

    script:
    // Cutoffs: calcuts' automatic choice, unless params.purge_dups_transition is set.
    // calcuts cannot place the haploid/diploid transition when the depth histogram
    // has no valley between the two peaks (one merged peak); it must then be given
    // (about 0.6x the depth peak). Lower/upper default to the automatic values.
    // With the transition unset the script is the original one, so existing runs stay cached.
    def calcuts_cmd = 'calcuts PB.stat > cutoffs'
    if (params.purge_dups_transition != null) {
        def lower = params.purge_dups_lower != null ? params.purge_dups_lower : '$(cut -f1 cutoffs_auto)'
        def upper = params.purge_dups_upper != null ? params.purge_dups_upper : '$(cut -f6 cutoffs_auto)'
        calcuts_cmd = [
            'calcuts PB.stat > cutoffs_auto 2> calcuts_auto.log',
            "calcuts -l ${lower} -m ${params.purge_dups_transition} -u ${upper} PB.stat > cutoffs",
        ].join('\n    ')
    }
    """
    ( set -o pipefail
      minimap2 -x map-ont -t ${task.cpus} ${assembly} ${reads} | gzip -c > reads_vs_asm.paf.gz )
    pbcstat reads_vs_asm.paf.gz
    ${calcuts_cmd}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purge_dups: "\$(purge_dups 2>&1 | grep -im1 version | sed 's/.*: *//' || echo unknown)"
        minimap2: "\$(minimap2 --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch PB.base.cov PB.stat cutoffs
    echo "PURGE_DUPS_COVERAGE stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purge_dups: "stub"
        minimap2: "stub"
    END_VERSIONS
    """
}

process PURGE_DUPS_RUN {
    label 'purge_dups_env'
    label 'res_purge_dups'
    tag   "purge_dups_run"

    conda 'bioconda::purge_dups=1.2.6 bioconda::minimap2'
    container 'quay.io/biocontainers/purge_dups:1.2.6--h577a1d6_3'
    publishDir "${params.outdir}/03b_purge_dups", mode: params.publish_mode

    input:
    path assembly
    path base_cov
    path cutoffs

    output:
    path "purged.fasta",    emit: purged_assembly
    path "haplotigs.fasta", emit: haplotigs
    path "dups.bed",        emit: dups_bed
    path "versions.yml", emit: versions

    script:
    """
    split_fa ${assembly} > split.fa
    ( set -o pipefail
      minimap2 -x asm5 -DP -t ${task.cpus} split.fa split.fa | gzip -c > split_self.paf.gz )

    purge_dups -2 -T ${cutoffs} -c ${base_cov} split_self.paf.gz > dups.bed

    get_seqs -e dups.bed ${assembly}
    # get_seqs writes purged.fa and hap.fa by default
    mv purged.fa purged.fasta
    mv hap.fa haplotigs.fasta
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purge_dups: "\$(purge_dups 2>&1 | grep -im1 version | sed 's/.*: *//' || echo unknown)"
        minimap2: "\$(minimap2 --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch purged.fasta haplotigs.fasta dups.bed
    echo "PURGE_DUPS_RUN stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purge_dups: "stub"
        minimap2: "stub"
    END_VERSIONS
    """
}

// purge_dups QC gate (bin/purge_qc_gate.py): BUSCO of the host bin before and after purging.
// FAIL stops the run (set the cutoffs by hand); WARN (duplication did not drop) is printed by main.nf.
process PURGE_QC_GATE {
    label 'res_light'
    tag   "purge_qc_gate"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/03b_purge_dups/qc", mode: params.publish_mode

    input:
    path before_summary, stageAs: 'before/*'
    path after_summary,  stageAs: 'after/*'
    path cutoffs

    output:
    path "purge_qc_report.tsv",   emit: report
    path "purge_qc_gate_mqc.tsv", emit: mqc
    path "purge_gate_status.txt", emit: status
    path "versions.yml",          emit: versions

    script:
    def source = params.purge_dups_transition != null ? 'params (purge_dups_lower/transition/upper)' : 'calcuts (automatic)'
    """
    purge_qc_gate.py --before ${before_summary} --after ${after_summary} --cutoffs ${cutoffs} \\
      --cutoffs-source '${source}' --busco-tolerance ${params.purge_busco_tolerance}

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version | sed 's/Python //')" > versions.yml
    """

    stub:
    """
    printf 'metric\\tbefore_purge\\tafter_purge\\tdelta\\tverdict\\n' > purge_qc_report.tsv
    cp purge_qc_report.tsv purge_qc_gate_mqc.tsv
    echo PASS > purge_gate_status.txt
    printf '"${task.process}":\\n    python: stub\\n' > versions.yml
    """
}

// Sequences shorter than params.min_contig_length (default 200 bp, the NCBI minimum) are removed;
// get_seqs -e leaves a few such fragments when it cuts contigs.
process LENGTH_FILTER {
    label 'res_light'
    tag   "length_filter"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/03b_purge_dups/length_filter", mode: params.publish_mode

    input:
    path assembly
    val  gate_status   // purge gate result (or 'skipped'): the filter runs only after the gate passed

    output:
    path "contigs.fasta",       emit: fasta
    path "length_removed.tsv",  emit: removed
    path "versions.yml",        emit: versions

    script:
    """
    length_filter.py ${assembly} contigs.fasta length_removed.tsv ${params.min_contig_length}
    test -s contigs.fasta

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version | sed 's/Python //')" > versions.yml
    """

    stub:
    """
    cp ${assembly} contigs.fasta
    printf 'name\\tlength\\n' > length_removed.tsv
    printf '"${task.process}":\\n    python: stub\\n' > versions.yml
    """
}
