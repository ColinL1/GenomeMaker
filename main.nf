// ============================================================
// main.nf — Coral Genome Assembly & Annotation Pipeline
// ============================================================
// params.busco_lineage — single lineage, e.g. 'anthozoa_odb12'.
// One BUSCO pass only, used for assembler comparison/ranking and for
// tagging host contigs during BlobToolKit filtering. No symbiont lineage.
// ============================================================

nextflow.enable.dsl = 2

include { PREPROCESS_READS }       from './subworkflows/00_preprocess_reads.nf'
include { KMER_GENOME_SIZE }       from './subworkflows/01_kmer_genome_size.nf'
include { MULTI_ASSEMBLER }        from './subworkflows/02_multi_assembler.nf'
include { SELECT_BEST_ASSEMBLY }   from './subworkflows/02b_select_best_assembly.nf'
include { PURGE_DUPS }             from './subworkflows/03b_purge_dups.nf'
include { ASSEMBLY_FILTERING }     from './subworkflows/03_assembly_filtering.nf'
include { ASSEMBLY_POLISHING }     from './subworkflows/04_assembly_polish.nf'
include { MITOGENOME_PIPELINE }    from './subworkflows/05_mitogenome.nf'
include { REPEATS_PIPELINE }       from './subworkflows/06_repeats.nf'
include { STRUCTURAL_ANNOTATION }  from './subworkflows/07_structural_annotation.nf'
include { FUNCTIONAL_ANNOTATION }  from './subworkflows/08_functional_annotation.nf'
include { GENERATE_SUMMARY; MULTIQC_FINAL } from './modules/09_summary.nf'
include { COPY_FILE }              from './modules/09.1_passthrough.nf'

workflow CORAL_GENOME_PIPELINE {

    def ont_reads_ch = channel.fromPath(params.ont_raw_reads, checkIfExists: true)

    def symbiont_path = params.symbiont_refs
        ? file(params.symbiont_refs.split(',').collect { s -> s.trim() }.find { s -> s }, checkIfExists: true)
        : null
    def symbiont_ch = channel.value(symbiont_path ?: [])

    def genomescope_path = (params.genomescope_script && file(params.genomescope_script).exists())
        ? file(params.genomescope_script)
        : null
    def genomescope_ch = channel.value(genomescope_path ?: [])

    def necat_config_path   = (params.necat_config && file(params.necat_config).exists()) ? file(params.necat_config) : null
    def necat_readlist_path = (params.necat_read_list && file(params.necat_read_list).exists()) ? file(params.necat_read_list) : null
    def necat_config_ch   = channel.value(necat_config_path ?: [])
    def necat_readlist_ch = channel.value(necat_readlist_path ?: [])

    def mito_ref_path = (params.mito_ref && file(params.mito_ref).exists()) ? file(params.mito_ref) : null

    def rnaseq_dir_path   = (params.rnaseq_dir && file(params.rnaseq_dir).exists()) ? file(params.rnaseq_dir) : null
    def protein_refs_path = (params.protein_refs && file(params.protein_refs).exists()) ? file(params.protein_refs) : null

    def funodb_path = (params.funannotate_db && file(params.funannotate_db).exists()) ? file(params.funannotate_db) : null
    def funodb_ch   = channel.value(funodb_path ?: [])

    // Optional MultiQC config (e.g. custom module order / titles). Falls back
    // to a harmless placeholder if params.multiqc_config isn't set, and the
    // MULTIQC process skips -c when it sees that placeholder name.
    def multiqc_config_path = (params.multiqc_config && file(params.multiqc_config).exists())
        ? file(params.multiqc_config)
        : file("NO_FILE")
    def multiqc_config_ch = channel.value(multiqc_config_path)

    // ===========================
    // 00. ONT READ PREPROCESSING
    // ===========================
    PREPROCESS_READS(ont_reads_ch, symbiont_ch)

    def porechopped_ch  = PREPROCESS_READS.out.reads_pc
    def assembly_reads  = PREPROCESS_READS.out.assembly_reads
    def polish_reads    = PREPROCESS_READS.out.polishing_reads

    // Read QC reports (NanoPlot/pycoQC/FastQC, etc.) emitted by the
    // preprocessing subworkflow.
    // ASSUMPTION: `00_preprocess_reads.nf` exposes a `qc_reports` output
    // channel with the raw read-QC report files. If your subworkflow uses
    // a different emit name, update the line below to match it; if it
    // doesn't expose one yet, add one there (e.g. emit the NanoPlot/
    // pycoQC/FastQC html+txt/json outputs) so MultiQC has read-level QC
    // to summarize. Until then, comment this line out and use
    // `channel.empty()` so the pipeline still runs.
    def read_qc_ch = PREPROCESS_READS.out.qc_dir_raw
        .mix(PREPROCESS_READS.out.qc_dir_trimmed)
        .flatMap { qc_dir -> file("${qc_dir}/**").findAll { it.isFile() } }
    // ===========================
    // 01. K-MER GENOME SIZE ESTIMATION
    // ===========================
    KMER_GENOME_SIZE(porechopped_ch, genomescope_ch)
    def genome_size_est = params.genome_size
        ? channel.value(params.genome_size)
        : KMER_GENOME_SIZE.out.genome_size

    // ===========================
    // 02. MULTI-ASSEMBLER DE NOVO ASSEMBLY
    // ===========================
    MULTI_ASSEMBLER(
        assembly_reads,
        necat_config_ch,
        necat_readlist_ch,
        genome_size_est,
        params.busco_lineage    // single lineage, e.g. anthozoa_odb12
    )

    def all_assemblies_ch      = MULTI_ASSEMBLER.out.assemblies
    def all_gfastats_ch        = MULTI_ASSEMBLER.out.gfastats
    def all_busco_ch           = MULTI_ASSEMBLER.out.busco_tables
    def busco_summary_ch       = MULTI_ASSEMBLER.out.busco_summary
    def comparison_ch          = MULTI_ASSEMBLER.out.comparison
    def comparison_json_ch     = MULTI_ASSEMBLER.out.comparison_json

    // ===========================
    // 02b. SELECT BEST ASSEMBLY (data-driven, not alphabetical)
    // Ranks by BUSCO completeness (penalizing duplication), N50 tie-break.
    // ===========================
    SELECT_BEST_ASSEMBLY(
        all_assemblies_ch,
        all_gfastats_ch,
        busco_summary_ch
    )

    def selected_assembly_path = SELECT_BEST_ASSEMBLY.out.selected_assembly
    def selected_name_ch       = SELECT_BEST_ASSEMBLY.out.selected_name
    def ranking_table_ch       = SELECT_BEST_ASSEMBLY.out.ranking_table

    selected_name_ch.view { name -> "Selected assembler: ${name}" }

    // ===========================
    // 02c. PURGE_DUPS
    // Remove duplicated haplotigs from the selected assembly using
    // ONT coverage, before contamination filtering / polishing.
    // ===========================
    PURGE_DUPS(selected_assembly_path, assembly_reads)
    def purged_assembly_path = PURGE_DUPS.out.purged_assembly
    def haplotigs_path       = PURGE_DUPS.out.haplotigs

    // ===========================
    // 03. ASSEMBLY FILTERING (BlobToolKit)
    // ===========================
    def meta_yaml_content =
    """
    sample: ${params.species_name}
    coverage: raw_ONT
    """
    def meta_yaml_file = file("${workDir}/meta_yaml/meta.yaml")
    meta_yaml_file.parent.mkdirs()
    meta_yaml_file.text = meta_yaml_content

    ASSEMBLY_FILTERING(
        purged_assembly_path,
        porechopped_ch,
        busco_summary_ch
            .combine(selected_name_ch)
            .filter { name, path, best -> name == best }
            .map    { name, path, best -> path }
            .first(),
        channel.value(meta_yaml_file)
    )

    def filtered_assembly_ch = ASSEMBLY_FILTERING.out.filtered_assembly
    def symbiont_assembly_ch = ASSEMBLY_FILTERING.out.symbiont_assembly

    // ===========================
    // 04. ASSEMBLY POLISHING
    // ===========================
    ASSEMBLY_POLISHING(
        filtered_assembly_ch,
        polish_reads,
        params.busco_lineage
    )

    def final_assembly_ch = ASSEMBLY_POLISHING.out.final_assembly
    def final_busco_ch    = ASSEMBLY_POLISHING.out.busco_tables

    // ===========================
    // 05. MITOGENOME EXTRACTION (optional)
    // ===========================
    def mito_bam_ch       = channel.empty()
    def mito_bai_ch       = channel.empty()
    def mito_coverage_ch  = channel.empty()
    def mito_flagstat_ch  = channel.empty()

    def run_mito = params.run_mitogenome && mito_ref_path
    if (run_mito) {
        MITOGENOME_PIPELINE(assembly_reads, channel.value(mito_ref_path))
        mito_bam_ch       = MITOGENOME_PIPELINE.out.mito_bam
        mito_bai_ch       = MITOGENOME_PIPELINE.out.mito_bai
        mito_coverage_ch  = MITOGENOME_PIPELINE.out.coverage
        mito_flagstat_ch  = MITOGENOME_PIPELINE.out.flagstat
    } else {
        log.info "Skipping mitogenome extraction (run_mitogenome=${params.run_mitogenome}, mito_ref=${params.mito_ref})"
    }

    // ===========================
    // 06. REPEAT ANNOTATION (optional)
    // ===========================
    def masked_assembly_ch = channel.empty()

    if (params.run_repeats) {
        REPEATS_PIPELINE(final_assembly_ch)
        masked_assembly_ch = REPEATS_PIPELINE.out.masked_assembly
    } else {
        COPY_FILE(final_assembly_ch)
        masked_assembly_ch = COPY_FILE.out.output
    }

    // ===========================
    // 07. STRUCTURAL ANNOTATION (optional)
    // ===========================
    def gff3_ch        = channel.empty()
    def proteins_faa_ch = channel.empty()

    def run_structural = params.run_structural_annotation && rnaseq_dir_path && protein_refs_path
    if (run_structural) {
        STRUCTURAL_ANNOTATION(masked_assembly_ch, rnaseq_dir_path, channel.value(protein_refs_path))
        gff3_ch        = STRUCTURAL_ANNOTATION.out.gff3
        proteins_faa_ch = STRUCTURAL_ANNOTATION.out.proteins
    } else {
        log.info "Skipping structural annotation (run_structural_annotation=${params.run_structural_annotation}, rnaseq_dir=${params.rnaseq_dir}, protein_refs=${params.protein_refs})"
    }

    // ===========================
    // 08. FUNCTIONAL ANNOTATION (optional)
    // ===========================
    def anno_dir_ch  = channel.empty()
    def stats_json_ch = channel.value([])

    def run_functional = params.run_functional_annotation && run_structural
    if (run_functional) {
        FUNCTIONAL_ANNOTATION(masked_assembly_ch, gff3_ch, proteins_faa_ch, funodb_ch)
        anno_dir_ch   = FUNCTIONAL_ANNOTATION.out.anno_dir
        stats_json_ch = FUNCTIONAL_ANNOTATION.out.stats_json
    } else {
        log.info "Skipping functional annotation (run_functional_annotation=${params.run_functional_annotation}, structural_annotation_ran=${run_structural})"
    }

    // ===========================
    // 09. SUMMARY REPORT
    // ===========================
    GENERATE_SUMMARY(
        all_gfastats_ch.map { name, p -> p },
        all_busco_ch.map { name, p -> p }.mix(final_busco_ch),
        comparison_ch,
        comparison_json_ch,
        stats_json_ch
    )

    // ===========================
    // 09b. MULTIQC — full-pipeline recap
    // Pulls together everything MultiQC can parse: pre-assembly read QC,
    // per-assembler gfastats/QUAST/BUSCO, the post-polish BUSCO re-run,
    // and (if generated) functional annotation stats. Optional steps that
    // didn't run simply contribute nothing, via channel.empty()/ifEmpty([]).
    // ===========================
    def multiqc_input_ch = channel.empty()
        .mix(read_qc_ch)
        .mix(all_gfastats_ch.map { name, p -> p })
        .mix(all_busco_ch.map { name, p -> p })
        .mix(final_busco_ch)
        .mix(comparison_ch)
        .mix(comparison_json_ch)
        .mix(stats_json_ch.filter { it })
        .flatten()
        .filter { it }
        .collect()
        .ifEmpty([])

    MULTIQC_FINAL(
        multiqc_input_ch,
        multiqc_config_ch
    )
}

workflow { CORAL_GENOME_PIPELINE() }
