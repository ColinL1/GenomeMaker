// ============================================================
// main.nf — Coral Genome Assembly & Annotation Pipeline
// ============================================================
// params.busco_lineage — single lineage, e.g. 'anthozoa_odb12'.
// One BUSCO pass only, used for assembler comparison/ranking and for
// tagging host contigs during BlobToolKit filtering. No symbiont lineage.
// ============================================================

nextflow.enable.dsl = 2

include { PREPROCESS_READS }       from './subworkflows/00_preprocess_reads.nf'
include { MERYL_READ_DB; KMER_GENOME_SIZE } from './subworkflows/01_kmer_genome_size.nf'
include { MULTI_ASSEMBLER }        from './subworkflows/02_multi_assembler.nf'
include { SELECT_BEST_ASSEMBLY }   from './subworkflows/02b_select_best_assembly.nf'
include { INPUT_ASSEMBLY_QC }      from './subworkflows/02c_input_assembly_qc.nf'
include { PURGE_DUPS }             from './subworkflows/03b_purge_dups.nf'
// Aliased so each checkpoint's outputs stay distinct.
include { ASSEMBLY_CHECKPOINT as CHECKPOINT_PURGED   } from './subworkflows/02d_assembly_checkpoint.nf'
include { ASSEMBLY_CHECKPOINT as CHECKPOINT_FILTERED } from './subworkflows/02d_assembly_checkpoint.nf'
include { ASSEMBLY_FILTERING }     from './subworkflows/03_assembly_filtering.nf'
include { ASSEMBLY_POLISHING }     from './subworkflows/04_assembly_polish.nf'
include { MITOGENOME_PIPELINE }    from './subworkflows/05_mitogenome.nf'
include { REPEATS_PIPELINE }       from './subworkflows/06_repeats.nf'
include { STRUCTURAL_ANNOTATION }  from './subworkflows/07_structural_annotation.nf'
include { FUNCTIONAL_ANNOTATION }  from './subworkflows/08_functional_annotation.nf'
include { RUN_SUMMARY }            from './subworkflows/09_summary.nf'
include { SCREEN_ASSEMBLY }        from './subworkflows/03c_screens.nf'
include { RAGTAG_SCAFFOLD }        from './modules/03d_scaffold.nf'
include { RELEASE }                from './subworkflows/10_release.nf'
include { linkage_evidence ; ncbi_linkage_terms } from './modules/10_release.nf'

workflow CORAL_GENOME_PIPELINE {

    // ---- Start-of-run checks: settings the late stages need, checked before any task runs
    def missing = []
    if (params.run_screens.toString().toBoolean()) {
        ['univec', 'fcs_gx_db', 'fcs_adaptor_image', 'fcs_gx_image', 'fcs_tax_id'].each { k ->
            if (!params[k]) missing << "--${k}"
        }
    }
    if (params.run_release.toString().toBoolean() && params.run_functional_annotation.toString().toBoolean()) {
        ['table2asn', 'ncbi_species', 'locus_tag', 'ncbi_sbt', 'ncbi_cmt', 'bioproject', 'biosample'].each { k -> if (!params[k]) missing << "--${k}" }
        // accessions go into every .sqn record: a wrong one attaches the genome to someone else's project
        if (params.bioproject && !(params.bioproject.toString() ==~ /PRJ[DEN][A-Z]\d+/)) {
            error "bioproject '${params.bioproject}' is not a BioProject accession (PRJNA..., PRJEB..., PRJDB...)"
        }
        if (params.biosample && !(params.biosample.toString() ==~ /SAM[DEN][A-Z]?\d+/)) {
            error "biosample '${params.biosample}' is not a BioSample accession (SAMN..., SAMEA..., SAMD...)"
        }
        if (params.uniprot_offline.toString().toBoolean() && !params.product_cache) {
            missing << '--product_cache (uniprot_offline takes product names from the cache only)'
        }
        if (params.ncbi_linkage_evidence && !(linkage_evidence() in ncbi_linkage_terms())) {
            error "ncbi_linkage_evidence '${params.ncbi_linkage_evidence}' is not an NCBI linkage-evidence term: " +
                  ncbi_linkage_terms().join(', ')
        }
        if (params.ncbi_species && params.ncbi_species.toString().contains('_')) {
            error "ncbi_species '${params.ncbi_species}' contains '_': use the organism name exactly as in the BioSample " +
                  "(NCBI cannot resolve it otherwise, and the record loses its lineage)"
        }
        if (!(params.ncbi_isolate || params.ncbi_strain)) {
            missing << '--ncbi_isolate or --ncbi_strain (NCBI requires one with a locus-tag prefix)'
        }
    }
    if (missing) {
        error "Missing settings: ${missing.join(', ')}.\n" +
              "Set them, or switch the stage off (--run_screens false / --run_release false)."
    }

    def ont_reads_ch = channel.fromPath(params.ont_raw_reads, checkIfExists: true)

    def symbiont_path = params.symbiont_refs
        ? file(params.symbiont_refs.split(',').collect { s -> s.trim() }.find { s -> s }, checkIfExists: true)
        : null
    def symbiont_ch = channel.value(symbiont_path ?: [])

    def necat_config_path   = (params.necat_config && file(params.necat_config).exists()) ? file(params.necat_config) : null
    def necat_readlist_path = (params.necat_read_list && file(params.necat_read_list).exists()) ? file(params.necat_read_list) : null
    def necat_config_ch   = channel.value(necat_config_path ?: [])
    def necat_readlist_ch = channel.value(necat_readlist_path ?: [])
    def necat_bin_path = (params.necat_path_bin && file(params.necat_path_bin).exists()) ? file(params.necat_path_bin) : null
    def necat_bin_ch   = channel.value(necat_bin_path ?: [])

    def mito_ref_path = (params.mito_ref && file(params.mito_ref).exists()) ? file(params.mito_ref) : null

    def rnaseq_dir_path   = (params.rnaseq_dir && file(params.rnaseq_dir).exists()) ? file(params.rnaseq_dir) : null
    def protein_refs_path = (params.protein_refs && file(params.protein_refs).exists()) ? file(params.protein_refs) : null

    def funodb_path = (params.funannotate_db && file(params.funannotate_db).exists()) ? file(params.funannotate_db) : null
    def funodb_ch   = channel.value(funodb_path ?: [])

    def eggnogdb_path = (params.eggnog_db && file(params.eggnog_db).exists()) ? file(params.eggnog_db) : null
    def eggnogdb_ch   = channel.value(eggnogdb_path ?: [])

    // Reference DBs / bundled helper scripts referenced by absolute host path
    // inside container script blocks. These MUST be declared as Nextflow
    // `path` process inputs (not interpolated straight from params) so that
    // apptainer/singularity/docker's automatic bind-mounting (autoMounts)
    // knows to expose the real host path inside the container -- otherwise
    // only the task work directory is bound and these paths are invisible
    // to the container even though they resolve fine outside of Nextflow.
    def blast_db_path      = file(params.blast_db)
    def blast_taxdump_path = file("${blast_db_path.getParent()}/new_taxdump")
    def diamond_db_path    = file(params.diamond_db)
    def trimmomatic_adapters_path = file(params.trimmomatic_adapters)

    // BUSCO offline lineage dir is optional; fall back to a placeholder
    // that Nextflow can stage even though it doesn't exist -- BUSCO_ASSEMBLY
    // only actually references it when params.busco_offline_dir is set.
    def busco_offline_path = (params.busco_offline_dir && file(params.busco_offline_dir).exists())
        ? file(params.busco_offline_dir)
        : file("${projectDir}/assets/NO_FILE", checkIfExists: true)
    def busco_offline_ch = channel.value(busco_offline_path)

    // Optional MultiQC config (e.g. custom module order / titles). Falls back
    // to a harmless placeholder if params.multiqc_config isn't set, and the
    // MULTIQC process skips -c when it sees that placeholder name.
    def multiqc_config_path = (params.multiqc_config && file(params.multiqc_config).exists())
        ? file(params.multiqc_config)
        : file("${projectDir}/assets/NO_FILE", checkIfExists: true)
    def multiqc_config_ch = channel.value(multiqc_config_path)

    // Fasta-assembly input mode: skip read preprocessing that only exists to
    // build a de novo assembly (k-mer genome size estimation, the
    // multi-assembler comparison, and best-assembly selection). ONT reads
    // are still required -- purge_dups coverage, BlobToolKit filtering,
    // polishing, and mitogenome extraction all need them.
    def use_input_assembly = params.input_assembly as boolean
    def input_assembly_path = use_input_assembly
        ? file(params.input_assembly, checkIfExists: true)
        : null

    // ===========================
    // 00. ONT READ PREPROCESSING
    // ===========================
    PREPROCESS_READS(ont_reads_ch, symbiont_ch)

    def porechopped_ch  = PREPROCESS_READS.out.reads_pc
    def assembly_reads  = PREPROCESS_READS.out.assembly_reads
    def polish_reads    = PREPROCESS_READS.out.polishing_reads
    // High-quality subset, k-mer genome sizing only (see CHOPPER_HQ_READS).
    def kmer_reads      = PREPROCESS_READS.out.kmer_reads

    // Raw, MultiQC-parseable read QC (FastQC zips, NanoPlot/NanoStat text) for
    // both the raw and trimmed stages.
    def read_qc_ch = PREPROCESS_READS.out.raw_qc_files

    // Every process reports its own tool version at runtime; those files are
    // accumulated here and collated in RUN_SUMMARY. Conditional stages mix inside
    // their branch, so a skipped stage simply contributes nothing.
    def versions_ch = PREPROCESS_READS.out.versions

    // ===========================
    // 01a. READ K-MER DATABASE
    // Built outside the genome-sizing branch because TWO consumers need it:
    // GenomeScope2 (de novo path only) and Merqury polishing QV (both paths,
    // including --input_assembly). Skipped only when neither is active.
    // ===========================
    def needs_sizing = !use_input_assembly && params.run_genomescope.toString().toBoolean() && !params.genome_size
    def needs_meryl  = needs_sizing || params.run_polish_qc.toString().toBoolean()

    def meryl_db_ch   = channel.value([])
    def meryl_hist_ch = channel.empty()
    if (needs_meryl) {
        MERYL_READ_DB(kmer_reads)
        meryl_db_ch   = MERYL_READ_DB.out.meryl_db
        meryl_hist_ch = MERYL_READ_DB.out.meryl_hist
        versions_ch   = versions_ch.mix(MERYL_READ_DB.out.versions)
    } else {
        log.info "Skipping read k-mer database (no genome sizing and run_polish_qc=false)"
    }

    // Channels the final MultiQC recap collects from optional/branching stages.
    // Declared up front so both pipeline paths can fill or leave them empty.
    def genomescope_out_ch = channel.empty()
    def quast_reports_ch   = channel.empty()
    def repeat_summary_ch  = channel.empty()
    def star_logs_ch       = channel.empty()
    def gff3_validation_ch = channel.empty()
    def compleasm_ch       = channel.empty()
    def omark_ch           = channel.empty()
    def ncrna_ch           = channel.empty()
    def stringtie_ch       = channel.empty()

    def all_gfastats_ch
    def all_busco_ch
    def busco_summary_ch
    def comparison_ch
    def comparison_json_ch
    def selected_assembly_path
    def selected_name_ch
    def ranking_table_ch = channel.empty()

    if (use_input_assembly) {
        log.info "Using supplied fasta assembly (input_assembly=${params.input_assembly}); skipping k-mer genome size estimation and de novo assembly/selection"

        INPUT_ASSEMBLY_QC(
            channel.value(input_assembly_path),
            params.busco_lineage,    // single lineage, e.g. anthozoa_odb12
            busco_offline_ch
        )

        all_gfastats_ch        = INPUT_ASSEMBLY_QC.out.gfastats
        all_busco_ch           = INPUT_ASSEMBLY_QC.out.busco_tables
        busco_summary_ch       = INPUT_ASSEMBLY_QC.out.busco_summary
        comparison_ch          = INPUT_ASSEMBLY_QC.out.comparison
        comparison_json_ch     = INPUT_ASSEMBLY_QC.out.comparison_json
        selected_assembly_path = INPUT_ASSEMBLY_QC.out.selected_assembly
        selected_name_ch       = INPUT_ASSEMBLY_QC.out.selected_name
        versions_ch            = versions_ch.mix(INPUT_ASSEMBLY_QC.out.versions)
    } else {
        // ===========================
        // 01. K-MER GENOME SIZE ESTIMATION
        // ===========================
        // Genome size: an explicit params.genome_size always wins; otherwise
        // estimate it with GenomeScope2 if run_genomescope is on. With neither,
        // the assemblers run without a size hint (Flye and hifiasm both accept
        // that) -- warn, because it changes their behaviour.
        def genome_size_est
        if (params.genome_size) {
            log.info "Using supplied genome_size=${params.genome_size}; skipping k-mer estimation"
            genome_size_est = channel.value(params.genome_size)
        }
        else if (params.run_genomescope.toString().toBoolean()) {
            KMER_GENOME_SIZE(meryl_hist_ch)
            genome_size_est    = KMER_GENOME_SIZE.out.genome_size
            genomescope_out_ch = KMER_GENOME_SIZE.out.genomescope_out
            versions_ch        = versions_ch.mix(KMER_GENOME_SIZE.out.versions)
        }
        else {
            log.warn "run_genomescope=false and no params.genome_size: assemblers will run without a genome-size hint"
            genome_size_est = channel.value([])
        }

        // ===========================
        // 02. MULTI-ASSEMBLER DE NOVO ASSEMBLY
        // ===========================
        MULTI_ASSEMBLER(
            assembly_reads,
            necat_config_ch,
            necat_readlist_ch,
            genome_size_est,
            params.busco_lineage,    // single lineage, e.g. anthozoa_odb12
            necat_bin_ch,
            busco_offline_ch
        )

        def all_assemblies_ch = MULTI_ASSEMBLER.out.assemblies
        all_gfastats_ch        = MULTI_ASSEMBLER.out.gfastats
        all_busco_ch           = MULTI_ASSEMBLER.out.busco_tables
        busco_summary_ch       = MULTI_ASSEMBLER.out.busco_summary
        comparison_ch          = MULTI_ASSEMBLER.out.comparison
        comparison_json_ch     = MULTI_ASSEMBLER.out.comparison_json
        quast_reports_ch       = MULTI_ASSEMBLER.out.quast_reports
        versions_ch            = versions_ch.mix(MULTI_ASSEMBLER.out.versions)

        // ===========================
        // 02b. SELECT BEST ASSEMBLY (data-driven, not alphabetical)
        // Ranks by BUSCO completeness (penalizing duplication), N50 tie-break.
        // ===========================
        SELECT_BEST_ASSEMBLY(
            all_assemblies_ch,
            all_gfastats_ch,
            busco_summary_ch,
            MULTI_ASSEMBLER.out.quast_reports,
            genome_size_est
        )

        selected_assembly_path = SELECT_BEST_ASSEMBLY.out.selected_assembly
        selected_name_ch       = SELECT_BEST_ASSEMBLY.out.selected_name
        ranking_table_ch       = SELECT_BEST_ASSEMBLY.out.ranking_table
        versions_ch            = versions_ch.mix(SELECT_BEST_ASSEMBLY.out.versions)
    }

    selected_name_ch.view { name -> "Selected assembler: ${name}" }

    // ===========================
    // 03b. PURGE_DUPS
    // Remove duplicated haplotigs from the selected assembly using
    // ONT coverage, before contamination filtering / polishing.
    // ===========================
    PURGE_DUPS(selected_assembly_path, assembly_reads)
    def purged_assembly_path = PURGE_DUPS.out.purged_assembly
    def haplotigs_path       = PURGE_DUPS.out.haplotigs
    versions_ch = versions_ch.mix(PURGE_DUPS.out.versions)

    // ===========================
    // 03b. CHECKPOINT -- measure the assembly AFTER purge_dups (02d_assembly_checkpoint)
    //
    // This serves two purposes at once, which is why it is not optional:
    //   1. QC. purge_dups deletes sequence; without a measurement here there is
    //      no way to tell a correct haplotig purge from an over-aggressive one.
    //   2. CORRECTNESS. BLOBTOOLS_ADD_BUSCO must be given a BUSCO table computed
    //      on the SAME assembly the BlobDir was built from: `get_seqs -e` renames
    //      split contigs, so a pre-purge table matches nothing -- and blobtools
    //      drops unmatched rows SILENTLY, leaving an empty BUSCO track.
    // ===========================
    CHECKPOINT_PURGED(
        purged_assembly_path,
        'purged',
        params.busco_lineage,
        busco_offline_ch
    )
    versions_ch = versions_ch.mix(CHECKPOINT_PURGED.out.versions)

    // ===========================
    // 03. ASSEMBLY FILTERING (BlobToolKit)
    // ===========================
    ASSEMBLY_FILTERING(
        purged_assembly_path,
        porechopped_ch,
        // The full_table.tsv (not short_summary.txt, which has no per-contig rows)
        // computed on the PURGED assembly -- the same sequences the BlobDir is built
        // from. Any other table gives an empty BUSCO track with no warning.
        CHECKPOINT_PURGED.out.busco_table.first(),
        blast_db_path,
        blast_taxdump_path,
        diamond_db_path
    )

    def filtered_assembly_ch = ASSEMBLY_FILTERING.out.filtered_assembly
    def symbiont_assembly_ch = ASSEMBLY_FILTERING.out.symbiont_assembly
    versions_ch = versions_ch.mix(ASSEMBLY_FILTERING.out.versions)

    // ===========================
    // 03. CHECKPOINT -- measure the assembly AFTER contamination filtering (02d_assembly_checkpoint)
    //
    // Purely QC, so it is behind a toggle. Its value is the DELTA against the
    // purged checkpoint above: that difference is exactly how much sequence the
    // taxon filter removed, which is the number the methods section needs and the
    // one that would reveal an over-aggressive filter.
    // ===========================
    def filtered_busco_ch    = channel.empty()
    def filtered_gfastats_ch = channel.empty()
    if (params.run_assembly_checkpoints.toString().toBoolean()) {
        CHECKPOINT_FILTERED(
            filtered_assembly_ch.first(),
            'filtered',
            params.busco_lineage,
            busco_offline_ch
        )
        filtered_busco_ch    = CHECKPOINT_FILTERED.out.busco_summary
        filtered_gfastats_ch = CHECKPOINT_FILTERED.out.gfastats
        versions_ch = versions_ch.mix(CHECKPOINT_FILTERED.out.versions)
    } else {
        log.info "Skipping the post-filtering QC checkpoint (run_assembly_checkpoints=false)"
    }

    // ===========================
    // 03c. CONTAMINATION SCREENS (FCS-adaptor, VecScreen, FCS-GX) on the host-filtered
    // contigs, before scaffolding. Unclear hits stop the run for review
    // (params.screen_decisions). The sequence is not changed after this step.
    // ===========================
    def screened_ch  = filtered_assembly_ch
    def screen_qc_ch = channel.empty()
    if (params.run_screens.toString().toBoolean()) {
        SCREEN_ASSEMBLY(filtered_assembly_ch.first())
        screened_ch  = SCREEN_ASSEMBLY.out.screened
        screen_qc_ch = SCREEN_ASSEMBLY.out.qc_files
        versions_ch  = versions_ch.mix(SCREEN_ASSEMBLY.out.versions)
    } else {
        log.warn "Contamination screens SKIPPED (run_screens=false); the release QC reports the screen checks as 'not run'"
    }

    // ===========================
    // 03d. SCAFFOLDING (RagTag), when params.scaffold_reference is set: after purge_dups
    // and the screens, so haplotigs and contaminants are never placed into chromosomes.
    // ===========================
    def scaffolds_ch = screened_ch
    if (params.scaffold_reference) {
        RAGTAG_SCAFFOLD(screened_ch, file(params.scaffold_reference, checkIfExists: true))
        scaffolds_ch = RAGTAG_SCAFFOLD.out.scaffolds
        versions_ch  = versions_ch.mix(RAGTAG_SCAFFOLD.out.versions)
    }

    // ===========================
    // 04. ASSEMBLY POLISHING
    // ===========================
    ASSEMBLY_POLISHING(
        scaffolds_ch,
        polish_reads,
        params.busco_lineage,
        busco_offline_ch,
        meryl_db_ch
    )

    def final_assembly_ch = ASSEMBLY_POLISHING.out.final_assembly
    def final_busco_ch    = ASSEMBLY_POLISHING.out.busco_tables
    versions_ch = versions_ch.mix(ASSEMBLY_POLISHING.out.versions)

    // ===========================
    // 05. MITOGENOME EXTRACTION (optional)
    // ===========================
    def mito_bam_ch       = channel.empty()
    def mito_bai_ch       = channel.empty()
    def mito_coverage_ch  = channel.empty()
    def mito_flagstat_ch  = channel.empty()
    def mito_fasta_ch     = channel.empty()
    def mito_gff_ch       = channel.empty()
    def mito_check_ch     = channel.empty()

    def run_mito = params.run_mitogenome.toString().toBoolean() && mito_ref_path
    if (run_mito) {
        // Optional Illumina shotgun DNA of the same individual (comma-separated
        // lists, R1 and R2 in the same order) -> GetOrganelle route.
        def mito_illumina_r1 = params.illumina_dna_r1 ? params.illumina_dna_r1.toString().split(',').collect { p -> file(p.trim(), checkIfExists: true) } : []
        def mito_illumina_r2 = params.illumina_dna_r2 ? params.illumina_dna_r2.toString().split(',').collect { p -> file(p.trim(), checkIfExists: true) } : []
        if (mito_illumina_r1.size() != mito_illumina_r2.size()) {
            error "illumina_dna_r1 and illumina_dna_r2 must list the same number of files, in the same order"
        }
        // `--mito_annot_ref_fasta false` on the command line arrives as the string 'false'
        def use_mito_annot     = params.mito_annot_ref_fasta && !(params.mito_annot_ref_fasta.toString() in ['false', 'null', ''])
        def mito_annot_ref     = use_mito_annot ? file(params.mito_annot_ref_fasta, checkIfExists: true) : []
        def mito_annot_ref_gff = use_mito_annot ? file(params.mito_annot_ref_gff,   checkIfExists: true) : []
        MITOGENOME_PIPELINE(assembly_reads, channel.value(mito_ref_path), selected_assembly_path,
                            polish_reads, mito_illumina_r1, mito_illumina_r2,
                            (mito_illumina_r1 + mito_illumina_r2).collect { f -> f.name }.join(', '),
                            mito_annot_ref, mito_annot_ref_gff)
        mito_bam_ch       = MITOGENOME_PIPELINE.out.mito_bam
        mito_bai_ch       = MITOGENOME_PIPELINE.out.mito_bai
        mito_coverage_ch  = MITOGENOME_PIPELINE.out.coverage
        mito_flagstat_ch  = MITOGENOME_PIPELINE.out.flagstat
        mito_fasta_ch     = MITOGENOME_PIPELINE.out.mitogenome
        mito_gff_ch       = MITOGENOME_PIPELINE.out.mito_gff
        mito_check_ch     = MITOGENOME_PIPELINE.out.mito_check
        versions_ch       = versions_ch.mix(MITOGENOME_PIPELINE.out.versions)
    } else {
        log.info "Skipping mitogenome extraction (run_mitogenome=${params.run_mitogenome}, mito_ref=${params.mito_ref})"
    }

    // ===========================
    // 06. REPEAT ANNOTATION (optional)
    // ===========================
    def masked_assembly_ch = channel.empty()
    def repeat_out_ch      = channel.empty()
    // Outputs with no MultiQC module, formatted as custom content by RUN_SUMMARY.
    def report_files_ch = ASSEMBLY_FILTERING.out.report_files
    if (run_mito) {
        report_files_ch = report_files_ch.mix(MITOGENOME_PIPELINE.out.report_files)
    }

    if (params.run_repeats.toString().toBoolean()) {
        REPEATS_PIPELINE(final_assembly_ch)
        masked_assembly_ch = REPEATS_PIPELINE.out.masked_assembly
        repeat_summary_ch  = REPEATS_PIPELINE.out.repeat_summary
        repeat_out_ch      = REPEATS_PIPELINE.out.repeat_tbl     // RepeatMasker .out
        report_files_ch    = report_files_ch.mix(REPEATS_PIPELINE.out.rrna_gff)
        versions_ch        = versions_ch.mix(REPEATS_PIPELINE.out.versions)
    } else {
        // masked_assembly_ch is consumed only by structural and functional
        // annotation, and both require run_repeats=true (enforced below), so with
        // repeats off it has no consumer.
        log.info "Skipping repeat annotation (run_repeats=false); structural/functional annotation unavailable"
    }

    // ===========================
    // 07. STRUCTURAL ANNOTATION (optional)
    // ===========================
    def gff3_ch        = channel.empty()
    def proteins_faa_ch = channel.empty()

    def run_structural = params.run_structural_annotation.toString().toBoolean() && rnaseq_dir_path && protein_refs_path

    // BRAKER is run with --softmasking, which requires a repeat-softmasked genome.
    // Repeats are therefore a PREREQUISITE for structural annotation: the
    // RepeatMasker -xsmall output is passed straight through, so masking is
    // computed once and reused rather than left to BRAKER.
    if (run_structural && !params.run_repeats.toString().toBoolean()) {
        error """
        run_structural_annotation=true requires run_repeats=true.

        BRAKER runs with --softmasking and needs the repeat-softmasked assembly
        that stage 06 produces. Without it BRAKER would receive an unmasked genome
        and predict genes across transposable elements.

        Either enable repeats (run_repeats=true) or disable structural annotation
        (run_structural_annotation=false).
        """.stripIndent()
    }

    if (run_structural) {
        STRUCTURAL_ANNOTATION(
            masked_assembly_ch,
            rnaseq_dir_path,
            channel.value(protein_refs_path),
            channel.value(trimmomatic_adapters_path),
            repeat_out_ch
        )
        gff3_ch            = STRUCTURAL_ANNOTATION.out.gff3
        proteins_faa_ch    = STRUCTURAL_ANNOTATION.out.proteins
        star_logs_ch       = STRUCTURAL_ANNOTATION.out.star_logs
        gff3_validation_ch = STRUCTURAL_ANNOTATION.out.gff3_validation
        compleasm_ch       = STRUCTURAL_ANNOTATION.out.compleasm
        omark_ch           = STRUCTURAL_ANNOTATION.out.omark
        ncrna_ch           = STRUCTURAL_ANNOTATION.out.ncrna
        stringtie_ch       = STRUCTURAL_ANNOTATION.out.stringtie
        report_files_ch    = report_files_ch.mix(
            STRUCTURAL_ANNOTATION.out.trnascan_out,
            STRUCTURAL_ANNOTATION.out.braker_results,
            STRUCTURAL_ANNOTATION.out.gff3,                 // refined.gff3
            STRUCTURAL_ANNOTATION.out.refined_validation,   // gt_refined.txt
            STRUCTURAL_ANNOTATION.out.trna_hc)
        versions_ch        = versions_ch.mix(STRUCTURAL_ANNOTATION.out.versions)
    } else {
        log.info "Skipping structural annotation (run_structural_annotation=${params.run_structural_annotation}, rnaseq_dir=${params.rnaseq_dir}, protein_refs=${params.protein_refs})"
    }

    // ===========================
    // 08. FUNCTIONAL ANNOTATION (optional)
    // ===========================
    def anno_dir_ch  = channel.empty()
    def stats_json_ch = channel.value([])

    def run_functional = params.run_functional_annotation.toString().toBoolean() && run_structural
    // The release needs the refined models and the functional evidence.
    def run_release    = params.run_release.toString().toBoolean() && run_functional
    if (run_functional) {
        FUNCTIONAL_ANNOTATION(masked_assembly_ch, gff3_ch, proteins_faa_ch, funodb_ch, eggnogdb_ch)
        anno_dir_ch   = FUNCTIONAL_ANNOTATION.out.anno_dir
        stats_json_ch = FUNCTIONAL_ANNOTATION.out.stats_json
        // Phobius may be a NO_ placeholder. With the release on, the report describes the
        // release's (NCBI) annotation instead of this first funannotate pass (below).
        report_files_ch = report_files_ch.mix(
            FUNCTIONAL_ANNOTATION.out.phobius.filter { f -> !f.name.startsWith('NO_') })
        if (!run_release) {
            report_files_ch = report_files_ch.mix(FUNCTIONAL_ANNOTATION.out.anno_dir)
        }
        versions_ch   = versions_ch.mix(FUNCTIONAL_ANNOTATION.out.versions)
    } else {
        log.info "Skipping functional annotation (run_functional_annotation=${params.run_functional_annotation}, structural_annotation_ran=${run_structural})"
    }

    report_files_ch = report_files_ch.mix(screen_qc_ch)   // screen actions and review (03c)

    // ===========================
    // 10. RELEASE: NCBI submission file, QC gate, release folder (before the summary,
    // so the MultiQC report describes the released annotation)
    // ===========================
    if (run_release) {
        RELEASE(
            masked_assembly_ch.first(),
            ASSEMBLY_POLISHING.out.presort_assembly.first(),
            gff3_ch,
            STRUCTURAL_ANNOTATION.out.refine_splits,
            STRUCTURAL_ANNOTATION.out.refine_tags,
            STRUCTURAL_ANNOTATION.out.refine_report,
            STRUCTURAL_ANNOTATION.out.trna_pseudogenes,
            STRUCTURAL_ANNOTATION.out.protein_identity,
            proteins_faa_ch,
            STRUCTURAL_ANNOTATION.out.refined_validation,
            anno_dir_ch,
            FUNCTIONAL_ANNOTATION.out.fixedproducts,
            FUNCTIONAL_ANNOTATION.out.eggnog,
            FUNCTIONAL_ANNOTATION.out.iprscan,
            FUNCTIONAL_ANNOTATION.out.phobius,
            funodb_ch,
            REPEATS_PIPELINE.out.repeat_tbl.mix(REPEATS_PIPELINE.out.repeat_summary,
                                                REPEATS_PIPELINE.out.repeat_gff, REPEATS_PIPELINE.out.repeat_lib),
            REPEATS_PIPELINE.out.repeat_lib,
            REPEATS_PIPELINE.out.classify_report,
            ASSEMBLY_FILTERING.out.symbiont_assembly,
            mito_fasta_ch,
            mito_gff_ch,
            mito_check_ch,
            ASSEMBLY_POLISHING.out.busco_summaries.first(),
            ASSEMBLY_POLISHING.out.final_qv,
            screen_qc_ch,
            versions_ch
        )
        stats_json_ch   = RELEASE.out.anno.map { d -> file("${d}/annotate_results/*.stats.json") }
                              .map { l -> l instanceof List ? l[0] : l }
        report_files_ch = report_files_ch.mix(RELEASE.out.anno, RELEASE.out.summary,
                                              RELEASE.out.ncbi.flatten().filter { f -> f.name.endsWith('.val') })
        versions_ch     = versions_ch.mix(RELEASE.out.versions)
        // An empty upstream channel would skip the gate and the folder without failing the run.
        RELEASE.out.release.ifEmpty { error "The release was not built: an input of the QC gate or the release folder was empty" }
    }

    // ===========================
    // 09. SUMMARY + INTEGRATED MULTIQC REPORT
    // Routed through a subworkflow like every other stage.
    // ===========================
    def multiqc_input_ch = channel.empty()
        // 00 -- read QC, raw and trimmed, so adapter/quality trimming is visible
        .mix(read_qc_ch)
        // 01 -- k-mer genome size model
        .mix(genomescope_out_ch)
        // 02 -- per-assembler QUAST (a first-class MultiQC module), gfastats, BUSCO
        .mix(quast_reports_ch.map    { name, p -> p })
        .mix(all_gfastats_ch.map     { name, p -> p })
        .mix(all_busco_ch.map        { name, p -> p })
        .mix(busco_summary_ch.map    { name, p -> p })
        // 02b -- which assembler was chosen and why
        .mix(ranking_table_ch)
        .mix(comparison_ch)
        .mix(comparison_json_ch)
        // 04 -- polishing: BUSCO both sides (+ the final product's if the gate fell
        // back to the unpolished assembly), Merqury QV, and the gate verdict
        .mix(final_busco_ch)
        .mix(ASSEMBLY_POLISHING.out.busco_stage_summaries)
        .mix(ASSEMBLY_POLISHING.out.merqury_pre.map  { stage, p -> p })
        .mix(ASSEMBLY_POLISHING.out.merqury_post.map { stage, p -> p })
        .mix(ASSEMBLY_POLISHING.out.polish_qc_report)
        .mix(ASSEMBLY_POLISHING.out.polish_gate_mqc)
        // 03b/03 -- what purge_dups and the taxon filter actually removed (checkpoints)
        .mix(CHECKPOINT_PURGED.out.busco_summary)
        .mix(CHECKPOINT_PURGED.out.gfastats)
        .mix(filtered_busco_ch)
        .mix(filtered_gfastats_ch)
        // 03 -- ONT coverage alignment rate against the purged assembly
        .mix(ASSEMBLY_FILTERING.out.coverage_flagstat)
        // 05 -- mitogenome coverage
        .mix(mito_coverage_ch)
        .mix(mito_flagstat_ch)
        // 06 -- repeat content
        .mix(repeat_summary_ch)
        // 07 -- RNA-seq alignment rates and annotation validation
        .mix(star_logs_ch)
        .mix(gff3_validation_ch)
        // 07b -- annotation completeness. OMArk is the INDEPENDENT metric: BRAKER's
        // own BUSCO is optimised against BUSCO via --busco_lineage, so it is circular.
        .mix(compleasm_ch)
        .mix(omark_ch)
        .mix(ncrna_ch)
        .mix(stringtie_ch)
        // 08 -- functional annotation yield
        .mix(stats_json_ch.filter { it })
        // all stages -- software provenance


    RUN_SUMMARY(
        all_gfastats_ch.map { name, p -> p },
        all_busco_ch.map { name, p -> p }.mix(final_busco_ch),
        // polishing-stage BUSCO summaries + the gate decision, so the summary can
        // say which assembly is the final product
        ASSEMBLY_POLISHING.out.busco_stage_summaries.mix(ASSEMBLY_POLISHING.out.polish_decision),
        comparison_ch,
        comparison_json_ch,
        stats_json_ch,
        multiqc_input_ch,
        report_files_ch,
        multiqc_config_ch,
        versions_ch
    )

}

workflow { CORAL_GENOME_PIPELINE() }
