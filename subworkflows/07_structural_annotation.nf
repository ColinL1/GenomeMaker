// ============================================================
// 07_structural_annotation.nf
// Structural gene annotation pipeline:
// tRNAscan-SE → RNA-seq trimming → STAR mapping → BRAKER4 → AGAT merge
// ============================================================

nextflow.enable.dsl = 2

include { COMPLEASM ; OMARK ; INFERNAL_CMSCAN ; STRINGTIE_ASSEMBLE } from '../modules/07b_annotation_qc.nf'
include { TRNASCAN_SE ; RNASEQ_TRIM ; STAR_GENOME_INDEX ; STAR_MAPPING ; MERGE_BAMS ; BRAKER4_RUN ; MERGE_ANNOTATIONS ; EXTRACT_PROTEINS ; GT_GFF3_VALIDATE } from '../modules/07_structural_annotation.nf'
include { TRNA_HC_FILTER ; HAPLOTIG_LONGEST ; HAPLOTIG_CLASSIFY ; GAP_SEGMENTS ; GAP_CLASSIFY ; REFINE_GFF ; GT_VALIDATE_REFINED ; GFFREAD_UNREFINED ; PROTEIN_IDENTITY } from '../modules/07c_annotation_refine.nf'
include { ISOFORM_PREP ; DIAMOND_ISOFORMS ; ISOFORM_BEDCOV ; ISOFORM_FILTER } from '../modules/07c_annotation_refine.nf'
include { DIAMOND_BLASTP as DIAMOND_ALLVSALL ; DIAMOND_BLASTP as DIAMOND_SEGMENTS } from '../modules/07c_annotation_refine.nf'

// ============================================================
// Workflow: STRUCTURAL_ANNOTATION
// tRNAscan-SE → RNA-seq → STAR → BRAKER4 → AGAT merge
//
// RNA-seq handling: `rnaseq_dir` (a plain directory path, not a
// channel) is used to build a channel of paired-end fastq tuples
// via `channel.fromFilePairs`, one item per sample. RNASEQ_TRIM and
// STAR_MAPPING are then forked/parallelized per-sample by Nextflow
// itself (no in-process for-loop over libraries). The per-sample
// STAR BAMs are collected and merged by a single MERGE_BAMS task.
// ============================================================
workflow STRUCTURAL_ANNOTATION {

    take:
        assembly
        rnaseq_dir       // plain directory path containing *_R1.fastq.gz / *_R2.fastq.gz
        protein_refs
        trimmomatic_adapters
        repeat_out       // RepeatMasker .out (repeat overlap filter for ncRNAs)

    main:
        def paired_reads_ch = channel.fromFilePairs("${rnaseq_dir}/*_R{1,2}.fastq.gz", checkIfExists: true)

        def versions_ch = channel.empty()

        TRNASCAN_SE(assembly)
        RNASEQ_TRIM(paired_reads_ch, trimmomatic_adapters)
        STAR_GENOME_INDEX(assembly)
        STAR_MAPPING(STAR_GENOME_INDEX.out.star_index_dir.first(), RNASEQ_TRIM.out.trimmed_reads)
        MERGE_BAMS(STAR_MAPPING.out.sample_bam.collect())

        // BRAKER4 takes the merged STAR BAM; it reads the transcript strand from the
        // XS tags, so no strand split is needed (and it has no --stranded option).
        BRAKER4_RUN(assembly, MERGE_BAMS.out.rnaseq_bam, protein_refs)
        def ncrna_tblout = channel.empty()
        if (params.rfam_cm && file(params.rfam_cm).exists()) {
            // Stage the DIRECTORY, not the .cm file: cmscan reads the cmpress
            // indices (.i1f/.i1i/.i1m/.i1p) that sit next to it, and staging the
            // single file would leave them behind.
            INFERNAL_CMSCAN(
                assembly,
                channel.value(file(params.rfam_cm).parent),
                file(params.rfam_cm).name,
                params.rfam_clanin ? file(params.rfam_clanin).name : ''
            )
            versions_ch = versions_ch.mix(INFERNAL_CMSCAN.out.versions)
            ncrna_tblout = INFERNAL_CMSCAN.out.tblout
        } else {
            log.warn "Infernal ncRNA SKIPPED (rfam_cm=${params.rfam_cm}) -- the final annotation " +
                     "will contain protein-coding genes and tRNAs only: no rRNA, snRNA, snoRNA " +
                     "or miRNA genes (these come from Rfam)."
        }

        // ncRNAs/rRNAs enter the annotation BEFORE functional annotation, so they reach
        // funannotate's GFF3, GenBank and NCBI TBL/SQN alike. Rfam family.txt types
        // each family (Gene vs Cis-reg ...); placeholders when Rfam is not available.
        def rfam_family_ch = (params.rfam_family && file(params.rfam_family).exists())
            ? channel.value(file(params.rfam_family))
            : channel.value(file("${projectDir}/assets/NO_RFAM_FAMILY", checkIfExists: true))
        MERGE_ANNOTATIONS(
            BRAKER4_RUN.out.braker_gtf,
            TRNASCAN_SE.out.trnascan_gff,
            ncrna_tblout.ifEmpty(file("${projectDir}/assets/NO_NCRNA_TBLOUT", checkIfExists: true)),
            rfam_family_ch
        )

        // Strict GFF3 spec validation in its own container (GenomeTools is not in
        // the AGAT image, which is why the old inline `gt gff3validator` could only fail).
        GT_GFF3_VALIDATE(MERGE_ANNOTATIONS.out.gff3)

        // GATE: everything downstream takes the merged GFF3 only once the validator
        // has passed. Previously the validator was a side branch, so protein
        // extraction and functional annotation started on a GFF3 that then failed
        // validation. The file itself is unchanged (same path), only its release
        // waits for the validation report.
        def validated_gff3 = MERGE_ANNOTATIONS.out.gff3
            .combine(GT_GFF3_VALIDATE.out.report)
            .map { gff, _report -> gff }
            .first()

        // ---- Refinement (07c), BEFORE protein extraction: every consumer below
        // (proteins, OMArk, StringTie, functional annotation, NCBI files) gets the
        // refined, canonical model set.
        TRNA_HC_FILTER(TRNASCAN_SE.out.trnascan_out, TRNASCAN_SE.out.trnascan_tbl)
        versions_ch = versions_ch.mix(TRNA_HC_FILTER.out.versions)

        def braker_res = BRAKER4_RUN.out.braker_results.first()
        def haplotig_tags = channel.value(file("${projectDir}/assets/NO_HAPLOTIG_TAGS", checkIfExists: true))
        if (params.run_haplotig_scan.toString().toBoolean()) {
            HAPLOTIG_LONGEST(braker_res)
            DIAMOND_ALLVSALL(HAPLOTIG_LONGEST.out.faa.map { f -> ['allvsall', f, f, '--more-sensitive -e 1e-10 -k 50'] })
            HAPLOTIG_CLASSIFY(HAPLOTIG_LONGEST.out.genes, DIAMOND_ALLVSALL.out.hits)
            haplotig_tags = HAPLOTIG_CLASSIFY.out.tags
            versions_ch = versions_ch.mix(DIAMOND_ALLVSALL.out.versions)
        }
        def gap_tags = channel.value(file("${projectDir}/assets/NO_GAP_TAGS", checkIfExists: true))
        if (params.related_proteomes) {
            GAP_SEGMENTS(braker_res, assembly, channel.value(file(params.related_proteomes, checkIfExists: true)))
            DIAMOND_SEGMENTS(GAP_SEGMENTS.out.faa.combine(GAP_SEGMENTS.out.refs)
                .map { q, r -> ['segments', q, r, '--more-sensitive -e 1e-5 -k 25'] })
            GAP_CLASSIFY(GAP_SEGMENTS.out.junctions, DIAMOND_SEGMENTS.out.hits)
            gap_tags = GAP_CLASSIFY.out.tags
        } else {
            log.warn "Gap-fusion scan SKIPPED (related_proteomes not set): gene models fused across " +
                     "assembly gaps are neither tagged nor split"
        }

        def transcript_edits = channel.value(file(params.transcript_edits ?: "${projectDir}/assets/NO_TRANSCRIPT_EDITS", checkIfExists: true))
        REFINE_GFF(braker_res, validated_gff3, assembly, TRNA_HC_FILTER.out.ids, repeat_out,
                   haplotig_tags, gap_tags, transcript_edits)

        // Isoform evidence filter: needs the related proteomes (homolog evidence)
        def models_gff3 = REFINE_GFF.out.gff
        def splits_ch = REFINE_GFF.out.splits
        if (params.related_proteomes && params.isoform_filter.toString().toBoolean()) {
            ISOFORM_PREP(REFINE_GFF.out.gff, assembly)
            DIAMOND_ISOFORMS(ISOFORM_PREP.out.faa, GAP_SEGMENTS.out.refs)
            ISOFORM_BEDCOV(ISOFORM_PREP.out.bed, MERGE_BAMS.out.rnaseq_bam, MERGE_BAMS.out.rnaseq_bai)
            ISOFORM_FILTER(REFINE_GFF.out.gff, assembly, braker_res, DIAMOND_ISOFORMS.out.hits,
                           ISOFORM_BEDCOV.out.cov, REFINE_GFF.out.splits)
            models_gff3 = ISOFORM_FILTER.out.gff
            splits_ch = ISOFORM_FILTER.out.splits
            versions_ch = versions_ch.mix(DIAMOND_ISOFORMS.out.versions)
        } else {
            log.warn "Isoform evidence filter SKIPPED (related_proteomes=${params.related_proteomes}, " +
                     "isoform_filter=${params.isoform_filter}): alternative isoforms are kept without an evidence check"
        }
        GT_VALIDATE_REFINED(models_gff3.mix(REFINE_GFF.out.pseudo))

        // GATE as above: the refined GFF3 is released downstream only after both
        // validations passed.
        def refined_gff3 = models_gff3
            .combine(GT_VALIDATE_REFINED.out.report.collect().map { _r -> 'validated' })
            .map { gff, _ok -> gff }
            .first()

        EXTRACT_PROTEINS(refined_gff3, assembly)

        // G2: proteins equal the unrefined ones except the split genes and removed transcripts
        GFFREAD_UNREFINED(validated_gff3, assembly)
        PROTEIN_IDENTITY(EXTRACT_PROTEINS.out.proteins_faa, GFFREAD_UNREFINED.out.proteins, splits_ch)

        // ---- Additive BRAKER4 steps: QC, ncRNA, transcript evidence ----
        // Each optional step is skipped when its database is absent, rather than
        // failing the run. See modules/07b_annotation_qc.nf for what braker.pl
        // already covers and is therefore NOT duplicated here.
        def compleasm_summary = channel.empty()
        if (params.run_polish_qc.toString().toBoolean()) {
            COMPLEASM(assembly, params.busco_lineage,
                      channel.value(params.compleasm_db
                          ? file(params.compleasm_db)
                          : file("${projectDir}/assets/NO_COMPLEASM_DB", checkIfExists: true)))
            versions_ch = versions_ch.mix(COMPLEASM.out.versions)
            compleasm_summary = COMPLEASM.out.summary
        }

        def omark_out = channel.empty()
        if (params.omamer_db && file(params.omamer_db).exists()) {
            OMARK(EXTRACT_PROTEINS.out.proteins_faa, refined_gff3, channel.value(file(params.omamer_db)))
            versions_ch = versions_ch.mix(OMARK.out.versions)
            omark_out = OMARK.out.outdir
        } else {
            log.warn "OMArk SKIPPED (omamer_db=${params.omamer_db}). This matters: BRAKER is run " +
                     "with --busco_lineage, which uses compleasm to pick the TSEBRA mode that " +
                     "MAXIMISES BUSCO -- so its BUSCO score is optimised against the metric and " +
                     "CANNOT be reported as an unbiased completeness figure. With OMArk off, " +
                     "nothing independently measures annotation completeness or contamination."
        }

        def stringtie_stats = channel.empty()
        if (params.run_stringtie.toString().toBoolean()) {
            STRINGTIE_ASSEMBLE(MERGE_BAMS.out.rnaseq_bam, refined_gff3)
            versions_ch = versions_ch.mix(STRINGTIE_ASSEMBLE.out.versions)
            stringtie_stats = STRINGTIE_ASSEMBLE.out.stats
        }

        versions_ch = versions_ch.mix(
            TRNASCAN_SE.out.versions, RNASEQ_TRIM.out.versions.first(),
            STAR_GENOME_INDEX.out.versions, STAR_MAPPING.out.versions.first(),
            MERGE_BAMS.out.versions,
            BRAKER4_RUN.out.versions, MERGE_ANNOTATIONS.out.versions, EXTRACT_PROTEINS.out.versions,
            GT_GFF3_VALIDATE.out.versions)

    emit:
        gff3     = refined_gff3                     // canonical (refined) models
        proteins = EXTRACT_PROTEINS.out.proteins_faa
        gff3_validation = GT_GFF3_VALIDATE.out.report
        refined_validation = GT_VALIDATE_REFINED.out.report
        trna_pseudogenes = REFINE_GFF.out.pseudo
        refine_splits    = splits_ch
        refine_tags      = REFINE_GFF.out.tags
        refine_report    = REFINE_GFF.out.report
        protein_identity = PROTEIN_IDENTITY.out.report
        trna_hc          = TRNA_HC_FILTER.out.hc_out
        star_logs       = STAR_MAPPING.out.star_log
        compleasm       = compleasm_summary
        omark           = omark_out
        ncrna           = ncrna_tblout
        stringtie       = stringtie_stats
        trnascan_out    = TRNASCAN_SE.out.trnascan_out
        braker_results  = BRAKER4_RUN.out.braker_results   // BRAKER4 QC (gene set, BUSCO, plots)
        versions        = versions_ch
}
