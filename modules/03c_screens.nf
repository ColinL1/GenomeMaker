// ============================================================
// 03c_screens.nf
// NCBI contamination screens and sequence hygiene:
//   FCS-adaptor, VecScreen (UniVec), FCS-GX -> review of unclear hits -> apply actions.
// Run on the host-filtered contigs (before scaffolding); FCS-adaptor and VecScreen
// also on the mitogenome (stage 05). Scripts: bin/vecscreen_classify.py,
// bin/screen_review.py, bin/apply_screen_actions.py.
// ============================================================

nextflow.enable.dsl = 2

// The FCS images ship without bash, which Nextflow needs inside a container, so they
// are run with `apptainer exec` from the host shell. Bind the task dir and the real
// locations of the inputs.
def fcs_binds(paths) {
    return "-B \$PWD " + paths.collect { p -> "-B \$(dirname \$(readlink -f ${p}))" }.join(' ')
}

process FCS_ADAPTOR {
    label 'res_mid_4h'
    tag   "${sname}"

    publishDir path: { "${params.outdir}/03c_screens/${sname}/fcs_adaptor" }, mode: params.publish_mode

    input:
    tuple val(sname), path(fasta)
    path image   // ncbi/fcs-adaptor .sif

    output:
    tuple val(sname), path("fcs_adaptor_report.txt"), emit: report
    path "fcs_adaptor.log"
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p out
    apptainer exec ${fcs_binds([fasta, image])} \$(readlink -f ${image}) \\
      /app/fcs/bin/av_screen_x -o \$PWD/out --euk \$(readlink -f ${fasta}) > run.log 2>&1
    cp out/fcs_adaptor_report.txt out/fcs_adaptor.log .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcs-adaptor: "${image.name}"
    END_VERSIONS
    """

    stub:
    """
    printf '#accession\\tlength\\taction\\trange\\tname\\n' > fcs_adaptor_report.txt
    touch fcs_adaptor.log
    printf '"${task.process}":\\n    fcs-adaptor: stub\\n' > versions.yml
    """
}

process VECSCREEN {
    tag   "${sname}"

    conda 'bioconda::blast=2.17.0'
    container 'quay.io/biocontainers/blast:2.17.0--h66d330f_0'

    input:
    tuple val(sname), path(fasta)
    path univec   // UniVec_Core FASTA

    output:
    tuple val(sname), path("raw.tsv"), emit: raw
    path "versions.yml", emit: versions

    script:
    // NCBI VecScreen parameters
    """
    makeblastdb -in ${univec} -dbtype nucl -out univec_core > /dev/null
    blastn -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true \\
      -evalue 700 -searchsp 1750000000000 -db univec_core -query ${fasta} -num_threads ${task.cpus} \\
      -outfmt "6 qseqid qstart qend sseqid sstart send pident length evalue bitscore score qlen" -out raw.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: "\$(blastn -version 2>&1 | head -1 | sed 's/blastn: //')"
        univec: "${univec.name}"
    END_VERSIONS
    """

    stub:
    """
    touch raw.tsv
    printf '"${task.process}":\\n    blast: stub\\n' > versions.yml
    """
}

process VECSCREEN_CLASSIFY {
    label 'res_light'
    tag   "${sname}"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir path: { "${params.outdir}/03c_screens/${sname}/vecscreen" }, mode: params.publish_mode

    input:
    tuple val(sname), path(raw)

    output:
    tuple val(sname), path("vecscreen_regions.tsv"), emit: regions
    path "vecscreen_summary.txt"

    script:
    """
    vecscreen_classify.py ${raw} vecscreen_regions.tsv > vecscreen_summary.txt
    """

    stub:
    """
    printf 'seq\\tstart\\tend\\tcategory\\tterminal\\tbest_subject\\tbest_score\\n' > vecscreen_regions.tsv
    echo stub > vecscreen_summary.txt
    """
}

process FCS_GX {
    tag "${sname}"

    publishDir path: { "${params.outdir}/03c_screens/${sname}/fcs_gx" }, mode: params.publish_mode

    input:
    tuple val(sname), path(fasta)
    path gxdb    // FCS-GX database directory (contains all.gxi, all.gxs, ...)
    path image   // ncbi/fcs-gx .sif

    output:
    tuple val(sname), path("*.fcs_gx_report.txt"), emit: report
    path "*.taxonomy.rpt"
    path "versions.yml", emit: versions

    script:
    """
    apptainer exec ${fcs_binds([fasta, image])} -B \$(readlink -f ${gxdb}) \\
      --env NCBI_FCS_REPORT_ANALYTICS=0 --env GX_NUM_CORES=${task.cpus} \$(readlink -f ${image}) \\
      run_gx --fasta \$(readlink -f ${fasta}) --out-dir \$PWD --gx-db \$(readlink -f ${gxdb})/all \\
      --tax-id ${params.taxon_id} > gx.log 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcs-gx: "${image.name}"
        fcs-gx-db: "${gxdb.name}"
    END_VERSIONS
    """

    stub:
    """
    printf '#seq_id\\tstart_pos\\tend_pos\\tseq_len\\taction\\tdiv\\tagg_cont_cov\\ttop_tax_name\\n' > ${sname}.${params.taxon_id}.fcs_gx_report.txt
    touch ${sname}.taxonomy.rpt
    printf '"${task.process}":\\n    fcs-gx: stub\\n' > versions.yml
    """
}

// Lists FCS-GX REVIEW rows and VecScreen Strong/Moderate matches, and turns the
// manual decisions (params.screen_decisions) into actions. This step always succeeds,
// so the list is published; APPLY_SCREEN_ACTIONS stops the run while any item is undecided.
process SCREEN_REVIEW {
    label 'res_light'
    tag   "${sname}"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir path: { "${params.outdir}/03c_screens/${sname}" }, mode: params.publish_mode

    input:
    tuple val(sname), path(gx_report), path(vec_regions)
    path decisions   // seq, start, end, decision, reason -- or a NO_ placeholder

    output:
    tuple val(sname), path("review/review_excluded.tsv"), path("review/review_mask.bed"), path("review/screen_review.tsv"), emit: actions
    tuple val(sname), path("review/accepted_fcs_gx.tsv"), path("review/accepted_vecscreen.tsv"), emit: accepted
    path "review/screen_review.tsv", emit: review

    script:
    def dec = decisions.name.startsWith('NO_') ? '' : "--decisions ${decisions}"
    """
    # exit status 3 = undecided items; reported by APPLY_SCREEN_ACTIONS after publishing
    screen_review.py --fcs-gx ${gx_report} --vecscreen ${vec_regions} ${dec} --outdir review || [ \$? -eq 3 ]
    """

    stub:
    """
    mkdir -p review
    touch review/review_excluded.tsv review/review_mask.bed review/screen_review.tsv
    printf 'seq\\tstart\\tend\\taction\\treason\\n' > review/accepted_fcs_gx.tsv
    printf 'seq\\tstart\\tend\\tcategory\\treason\\n' > review/accepted_vecscreen.tsv
    """
}

// Applies FCS actions (EXCLUDE/TRIM/FIX), the review decisions, and NCBI sequence
// hygiene (terminal N runs removed, N runs upper-cased).
process APPLY_SCREEN_ACTIONS {
    label 'res_mid_4h'
    tag   "${sname}"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir path: { "${params.outdir}/03c_screens/${sname}" }, mode: params.publish_mode,
        saveAs: { fn -> fn.endsWith('.fasta') ? "screened/${fn}" : fn }

    input:
    tuple val(sname), path(fasta), path(reports), path(review_excluded), path(mask_bed), path(review)

    output:
    path "${sname}.screened.fasta", emit: fasta
    path "screen_actions.tsv",      emit: log

    script:
    """
    n=\$(awk -F'\\t' 'NR>1 && \$7=="UNDECIDED"' ${review} | wc -l)
    if [ "\$n" -gt 0 ]; then
        echo "ERROR: \$n contamination-screen hits need a decision. See ${params.outdir}/03c_screens/${sname}/review/screen_review.tsv;" >&2
        echo "copy the UNDECIDED rows into the file given by --screen_decisions (columns seq, start, end," >&2
        echo "decision = exclude|mask|trim|keep, reason) and resume the run with -resume." >&2
        exit 1
    fi
    rep=""
    for r in ${reports} ${review_excluded}; do rep="\$rep --report \$r"; done
    mb=""
    if [ -s ${mask_bed} ]; then mb="--mask-bed ${mask_bed}"; fi
    apply_screen_actions.py --fasta ${fasta} \$rep \$mb \\
      --out-fasta ${sname}.screened.fasta --out-gff-dir lifted --out-log screen_actions.tsv
    test -s ${sname}.screened.fasta
    """

    stub:
    """
    cp ${fasta} ${sname}.screened.fasta
    touch screen_actions.tsv
    """
}
