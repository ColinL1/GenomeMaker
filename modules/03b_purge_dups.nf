// ============================================================
// 03b_purge_dups.nf
// Standard purge_dups pipeline (read-depth + self-alignment based
// haplotig purging), using ONT coverage. Runs on the SELECTED
// assembly, before BlobToolKit filtering.
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
