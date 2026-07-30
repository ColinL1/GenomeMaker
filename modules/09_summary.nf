// ============================================================
// 09_summary.nf
// Aggregate all QC results into a final summary report
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: GENERATE_SUMMARY
// Collect gfastats, BUSCO, QUAST, annotation stats into reports
// ============================================================
process GENERATE_SUMMARY {
    tag   "generate_summary"
    cpus  { 2 * task.attempt }
    memory { 8.GB * task.attempt }
    time   { 2.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'conda-forge::pandas conda-forge::jinja2'
    container 'quay.io/biocontainers/mulled-v2-5ff8b00c2d7f6173e034c115dfe295627ff99689:beb0ad5f49ec2904f79edffacd41bba38492e881-0'

    publishDir "${params.outdir}/summary", mode: params.publish_mode

    input:
    path gfastats_files
    path busco_files
    path comparison_tsv
    path comparison_json
    path annotation_stats

    output:
    path "summary/assembly_comparison.tsv", emit: summary_tsv
    path "summary/assembly_comparison.json", emit: summary_json
    path "summary/pipeline_report.txt",      emit: pipeline_report

    script:
    """
    mkdir -p summary

    # Copy comparison reports
    cp ${comparison_tsv} summary/assembly_comparison.tsv 2>/dev/null || true
    cp ${comparison_json} summary/assembly_comparison.json 2>/dev/null || true

    # Generate pipeline report
    cat > summary/pipeline_report.txt << 'REPORT'
    ============================================================
    Coral Genome Assembly & Annotation Pipeline — Summary Report
    ============================================================

    Assembly QC Summary
    -------------------
    (See assembly_comparison.tsv for detailed per-assembler stats)

    BUSCO Summary
    -------------
    (See assembly_comparison.tsv for BUSCO scores)

    Annotation Summary
    ------------------
    (See annotation_stats.json for details)

    ============================================================
    REPORT

    echo "Summary report generated"
    """

    stub:
    """
    mkdir -p summary
    touch summary/assembly_comparison.tsv summary/assembly_comparison.json summary/pipeline_report.txt
    echo "GENERATE_SUMMARY stub completed"
    """
}

// ============================================================
// Process: MULTIQC
// Aggregate every QC-producing step of the pipeline into a
// single MultiQC report: read QC (pre-assembly), per-assembler
// gfastats/QUAST/BUSCO, the post-polish BUSCO re-run, and
// functional annotation stats (e.g. funannotate/InterProScan/
// eggNOG summary json, if MultiQC-recognisable).
//
// This process is intentionally permissive about its inputs:
// it just stages whatever files it is given into one directory
// and lets `multiqc` auto-detect what it recognises. Missing or
// empty optional channels (e.g. mitogenome/annotation not run)
// are fine — MultiQC will simply report on what's present.
// ============================================================
process MULTIQC_FINAL {
    tag   "multiqc"
    cpus  { 2 * task.attempt }
    memory { 4.GB * task.attempt }
    time   { 1.h * task.attempt }
    errorStrategy 'retry'
    maxRetries params.max_retries

    conda 'bioconda::multiqc'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'

    publishDir "${params.outdir}/summary/multiqc", mode: params.publish_mode

    input:
    path qc_files, stageAs: 'qc_inputs/*'
    path multiqc_config

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data
    path "multiqc_report_data.zip", emit: data_zip, optional: true

    script:
    def config_arg = (multiqc_config && multiqc_config.name != 'NO_FILE') ? "-c ${multiqc_config}" : ""
    """
    multiqc qc_inputs \\
        --force \\
        --filename multiqc_report.html \\
        --outdir . \\
        ${config_arg}

    # Keep a zipped copy of the parsed data alongside the html report
    if [ -d multiqc_data ]; then
        zip -qr multiqc_report_data.zip multiqc_data
    fi
    """

    stub:
    """
    touch multiqc_report.html
    mkdir -p multiqc_data
    echo stub > multiqc_data/stub.txt
    touch multiqc_report_data.zip
    echo "MULTIQC stub completed"
    """
}
