// ============================================================
// 09_summary.nf
// Final reporting: software provenance + the integrated MultiQC report.
//
// Exists so main.nf goes through a subworkflow for EVERY stage. Stage 09 was the
// one exception -- main.nf imported modules/09_summary.nf and 09b directly, which
// broke the otherwise-uniform pattern and put report-assembly logic in main.nf.
// ============================================================

nextflow.enable.dsl = 2

include { GENERATE_SUMMARY; MULTIQC_CUSTOM_CONTENT; MULTIQC_FINAL } from '../modules/09_summary.nf'
include { COLLATE_VERSIONS }               from '../modules/09b_software_versions.nf'

workflow RUN_SUMMARY {

    take:
        gfastats_files      // collected gfastats reports
        busco_files         // collected BUSCO tables
        polish_files        // polishing-stage BUSCO summaries + polish_decision.txt
        comparison_tsv
        comparison_json
        annotation_stats
        multiqc_files       // everything MultiQC can parse, from every stage
        report_files        // outputs MultiQC has no module for (formatted as custom content)
        multiqc_config
        versions_files      // per-process versions.yml from every stage

    main:
        // collect(): these are per-file channels; without it GENERATE_SUMMARY ran once
        // per file pairing and each report saw only one gfastats/BUSCO file.
        GENERATE_SUMMARY(
            gfastats_files.collect().ifEmpty([]),
            busco_files.collect().ifEmpty([]),
            polish_files.collect().ifEmpty([]),
            comparison_tsv,
            comparison_json,
            annotation_stats
        )

        // Real per-tool versions, collated from what each process reported at runtime.
        // collectFile concatenates into ONE file. Staging the 63 per-process files
        // individually fails with "input file name collision" -- they are all
        // literally named versions.yml, and stageAs preserves basenames.
        // sort by CONTENT, not by path: several stages emit `.versions.first()`, so
        // which task's versions.yml arrives (and hence its path) differs between runs
        // and a path sort reordered the file on every -resume.
        COLLATE_VERSIONS(
            versions_files.collectFile(name: 'collated_versions.yml', sort: { f -> f.text })
        )

        // The versions table is part of the report, so it joins the input set here
        // rather than being mixed in by the caller.
        // Tools without a MultiQC module -> *_mqc custom content (bin/multiqc_custom_content.py).
        MULTIQC_CUSTOM_CONTENT(
            multiqc_files.mix(report_files).flatten().filter { it }.unique().collect().ifEmpty([])
        )

        MULTIQC_FINAL(
            multiqc_files.mix(COLLATE_VERSIONS.out.mqc, MULTIQC_CUSTOM_CONTENT.out.mqc).flatten().filter { it }.collect().ifEmpty([]),
            multiqc_config
        )

    emit:
        report            = MULTIQC_FINAL.out.report
        custom_content    = MULTIQC_CUSTOM_CONTENT.out.mqc   // *_mqc tables (the QC gate reads the funannotate one)
        multiqc_data      = MULTIQC_FINAL.out.data
        software_versions = COLLATE_VERSIONS.out.table
        pipeline_report   = GENERATE_SUMMARY.out.pipeline_report
}
