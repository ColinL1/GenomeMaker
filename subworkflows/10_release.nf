// ============================================================
// 10_release.nf (subworkflow)
// The NCBI-ready release, from the finished assembly and its one annotation run:
// product curation -> funannotate NCBI pass -> .tbl fixes -> table2asn preflight -> what table2asn
// changed (.ecn, .fixedproducts) fed back into the curation and the .tbl -> table2asn (+ mitogenome),
// scaffold name map, release files synced to the submitted .tbl, high-confidence subset, symbiont BUSCO,
// summary + README, QC gate (a FAIL stops the run), release folder <outdir>/release/<release_name>/.
// The released sequence is the final masked assembly as it is.
// ============================================================

nextflow.enable.dsl = 2

include { CURATE_PRODUCTS ; CURATE_PRODUCTS as CURATE_PRODUCTS_FINAL ; MAKE_FIX_FILES ; FUNANNOTATE_NCBI ;
          FIX_TBL as FIX_TBL_PREFLIGHT ; FIX_TBL ; MITO_TBL ; TABLE2ASN as TABLE2ASN_PREFLIGHT ; TABLE2ASN ;
          SCAFFOLD_NAME_MAP ; RELEASE_SETS ; SYMBIONT_BUSCO ; RELEASE_METHODS ; RELEASE_REPORT ; RELEASE_MQC_TABLE ;
          ANNOTATION_QC_GATE ; RELEASE_ASSEMBLE ; linkage_evidence } from '../modules/10_release.nf'

// Placeholder for an absent optional input (assets/NO_*)
def no_file(name) {
    return file("${projectDir}/assets/${name}", checkIfExists: true)
}

// Optional file param as a value channel, or its NO_ placeholder
def opt_file(value, placeholder) {
    return channel.value(value ? file(value, checkIfExists: true) : no_file(placeholder))
}

// Files of one extension from a process's output list, or the NO_ placeholder when there are none
def by_ext(ch, ext, placeholder) {
    return ch.flatten().filter { f -> f.name.endsWith(ext) }.collect().ifEmpty([no_file(placeholder)])
}

// Methods section of the release README, from the steps and settings this run used. {{tool}} is filled
// with the version the run's processes recorded for that tool (RELEASE_METHODS, bin/release_methods.py).
def release_methods() {
    def purge = params.purge_dups_transition != null
        ? "calcuts transition ${params.purge_dups_transition} (lower ${params.purge_dups_lower ?: 'automatic'}, upper ${params.purge_dups_upper ?: 'automatic'})"
        : 'automatic calcuts cutoffs'
    def before_scaffolding = params.scaffold_reference ? ' on the contigs, before scaffolding' : ''
    def lines = ["- Haplotig purging: purge_dups {{purge_dups}}${before_scaffolding}, ${purge}; `get_seqs -e`."]
    if (!params.skip_blobtools.toString().toBoolean()) {
        lines << "- Host filter: BlobToolKit {{blobtoolkit}}${before_scaffolding}; the symbiont bin is set aside."
    }
    if (params.run_screens.toString().toBoolean()) {
        lines << "- Contamination: FCS-adaptor {{fcs-adaptor}} (--euk), VecScreen (UniVec_Core, BLAST {{blast}}) and FCS-GX {{fcs-gx}} " +
                 "(tax-id ${params.fcs_tax_id})${before_scaffolding}; clear-cut actions applied, unclear hits reviewed " +
                 "(`qc/screen_review.tsv`). The sequence is not changed after this step. The mitogenome: FCS-adaptor and VecScreen."
    }
    if (params.scaffold_reference) {
        lines << "- Scaffolding: RagTag {{ragtag}} `scaffold` against `${file(params.scaffold_reference).name}`" +
                 (params.ragtag_infer_gaps.toString().toBoolean() ? ' (gap sizes inferred)' : ' (100 bp gaps of unknown size)') + '.'
    }
    lines += [
        "- Repeats: RepeatModeler {{repeatmodeler}} + EDTA {{edta}} library; Unknown families classified with RepeatClassifier; " +
            "RepeatMasker {{repeatmasker}} soft-masking with the classified library.",
        "- Gene models: BRAKER4 {{braker4}} (RNA-seq + proteins), tRNAscan-SE {{trnascan-se}} high-confidence set, " +
            "Rfam/Infernal {{infernal}} ncRNAs; UTRs from BRAKER's stringtie2utr, gene membership by a frame-aware " +
            "CDS-overlap rule (>= ${params.refine_frame_min}), Rfam filters, evidence and haplotig tags" +
            (params.related_proteomes ? ", models fused across assembly gaps split where >= ${params.gap_split_min_species} related proteomes support it" : '') +
            (params.related_proteomes && params.isoform_filter.toString().toBoolean()
                ? "; alternative isoforms kept only with evidence of their own (an intron with >= ${params.isoform_min_reads} RNA-seq reads, " +
                  "RNA-seq depth, or better cover of a related protein), one isoform kept per gene" : '') + '.',
        "- Functional annotation: eggNOG-mapper {{eggnog-mapper}}" +
            (params.interproscan_data || params.interproscan_dir ? ', InterProScan {{interproscan}}' : '') +
            (params.run_phobius.toString().toBoolean() ? ', Phobius {{phobius}}' : '') + ", funannotate {{funannotate}}; product " +
            "names checked against UniProt reviewed entries by gene symbol (UniProt release {{uniprot}}; `product_curation.tsv`).",
        "- NCBI: table2asn {{table2asn}} with the organism '${params.ncbi_species}'" +
            (params.locus_tag ? ", locus-tag prefix ${params.locus_tag}" : '') + "; gaps U / 100 bp / ${linkage_evidence()}; " +
            "EC numbers and product names table2asn changed in a preflight run were applied to the submitted annotation.",
    ]
    return lines.join('\n')
}

workflow RELEASE {

    take:
        assembly            // final soft-masked assembly: annotated as it is, released as it is
        presort_assembly    // the same sequences before FUNANNOTATE_SORT renamed them (name map)
        refined_gff
        refine_splits
        refine_tags
        refine_report
        protein_identity
        refined_proteins
        gt_reports          // GT validation of the refined GFF3 files
        anno_pass1          // funannotate annotate (FUNCTIONAL_ANNOTATION)
        fixedproducts       // its *.fixedproducts (may be empty)
        eggnog
        iprscan             // or NO_ placeholder
        phobius             // or NO_ placeholder
        funannotate_db
        repeat_files        // RepeatMasker .out, .tbl, repeats.gff3, classified library
        repeat_lib          // classified library (G8)
        classify_report     // RepeatClassifier family counts
        symbiont_assembly   // may be empty (skip_blobtools, or no symbiont fraction)
        mito_fasta          // channel, may be empty
        mito_gff            // channel, may be empty
        mito_check          // channel, may be empty
        busco_summary       // BUSCO short summary of the final assembly (stage 04)
        qv                  // Merqury QV of the final assembly (stage 04), may be empty
        screen_files        // screen_actions.tsv + screen_review.tsv (03c), may be empty
        versions            // versions.yml of the stages before the release (README Methods)

    main:
        def versions_ch = channel.empty()
        def sbt = channel.value(file(params.ncbi_sbt, checkIfExists: true))
        def t2a = file(params.table2asn, checkIfExists: true)
        def cmt = channel.value(file(params.ncbi_cmt, checkIfExists: true))
        def draft = opt_file(params.product_draft, 'NO_PRODUCT_DRAFT')
        // passed as inputs: a param read inside a script block does not enter the task hash (-resume would reuse the old result)
        def lookup_all = params.product_lookup_all.toString().toBoolean() ? '--lookup-all' : ''
        def strip = params.strip_gene_suffix.toString().toBoolean() ? '--strip-gene-suffix' : ''
        def pass1_fixed = fixedproducts.collect().ifEmpty([no_file('NO_FIXEDPRODUCTS')])
        def prev_cache = file("${params.outdir}/10_release/product_curation/uniprot_cache.json")
        if (!params.product_cache && prev_cache.exists()) {
            log.info "Product curation queries UniProt again; to reuse the earlier lookups (same names), set --product_cache ${prev_cache}"
        }

        // ---- products, funannotate NCBI pass
        CURATE_PRODUCTS('pass1', lookup_all, anno_pass1, pass1_fixed, draft, opt_file(params.product_cache, 'NO_PRODUCT_CACHE'), funannotate_db)
        MAKE_FIX_FILES(refined_gff, assembly, CURATE_PRODUCTS.out.curation)
        FUNANNOTATE_NCBI(refined_gff, assembly, eggnog, iprscan, phobius, funannotate_db,
                         MAKE_FIX_FILES.out.fix, MAKE_FIX_FILES.out.remove, sbt)
        def anno_ch = FUNANNOTATE_NCBI.out.dir

        // ---- mitogenome (when Liftoff annotated it), in the same submission
        MITO_TBL(mito_fasta, mito_gff, mito_check)
        def mito_fsa = MITO_TBL.out.fsa.ifEmpty(no_file('NO_MITO_FSA')).first()
        def mito_tbl = MITO_TBL.out.tbl.ifEmpty(no_file('NO_MITO_TBL')).first()

        // ---- preflight: .tbl fixes, table2asn; what table2asn changed goes back into the curation and the .tbl
        FIX_TBL_PREFLIGHT('preflight', strip, anno_ch, refined_gff, no_file('NO_ECN'), no_file('NO_PRODUCT_CURATION'),
                          no_file('NO_FIXEDPRODUCTS'))
        TABLE2ASN_PREFLIGHT('preflight', FIX_TBL_PREFLIGHT.out.dir, t2a, sbt, cmt, mito_fsa, mito_tbl)
        def pre_ecn   = by_ext(TABLE2ASN_PREFLIGHT.out.files, '.ecn', 'NO_ECN')
        def pre_fixed = by_ext(TABLE2ASN_PREFLIGHT.out.files, '.fixedproducts', 'NO_FIXEDPRODUCTS')
        CURATE_PRODUCTS_FINAL('final', lookup_all, anno_pass1, fixedproducts.mix(pre_fixed).flatten().filter { f -> !f.name.startsWith('NO_') }
                                  .collect().ifEmpty([no_file('NO_FIXEDPRODUCTS')]),
                              draft, CURATE_PRODUCTS.out.cache, funannotate_db)

        // ---- the submission
        FIX_TBL('final', strip, anno_ch, refined_gff, pre_ecn, CURATE_PRODUCTS_FINAL.out.curation, pre_fixed)
        TABLE2ASN('final', FIX_TBL.out.dir, t2a, sbt, cmt, mito_fsa, mito_tbl)

        // ---- names, sets (synced to the submitted .tbl), symbiont BUSCO
        SCAFFOLD_NAME_MAP(presort_assembly, assembly)
        RELEASE_SETS(anno_ch, FIX_TBL.out.dir, TABLE2ASN.out.files.flatten().filter { f -> f.name.endsWith('.gbf') }.first(),
                     refined_gff, refine_tags, refine_splits, FIX_TBL.out.rrna)
        def sets        = RELEASE_SETS.out.files.flatten()
        def rel_gff     = sets.filter { f -> f.name.endsWith('.gff3') && !f.name.contains('.hc.') && f.name != 'trna_pseudogenes.gff3' }.first()
        def sets_report = sets.filter { f -> f.name == 'release_sets_report.tsv' }.first()
        def sym_summary = channel.value(no_file('NO_SYMBIONT_BUSCO'))
        if (params.symbiont_busco_lineage) {
            SYMBIONT_BUSCO(symbiont_assembly, file(params.symbiont_busco_lineage, checkIfExists: true))
            // no symbiont bin (skip_blobtools, or nothing assigned to the symbiont): no BUSCO, not a missing input
            sym_summary = SYMBIONT_BUSCO.out.summary.flatten().filter { f -> f.name.endsWith('.txt') }
                .ifEmpty {
                    log.info "Symbiont BUSCO not run: no symbiont bin"
                    no_file('NO_SYMBIONT_BUSCO')
                }.first()
            versions_ch = versions_ch.mix(SYMBIONT_BUSCO.out.versions)
        }
        versions_ch = versions_ch.mix(CURATE_PRODUCTS_FINAL.out.versions, FUNANNOTATE_NCBI.out.versions, TABLE2ASN.out.versions)

        // ---- summary + README, the MultiQC protein-coding table, QC gate, folder
        def template = channel.of(release_methods()).collectFile(name: 'methods_template.md', newLine: true)
        RELEASE_METHODS(template, versions.mix(versions_ch).collectFile(name: 'all_versions.yml', newLine: true))
        def repeat_tbl = repeat_files.flatten().filter { f -> f.name.endsWith('.tbl') }.first()
        def qv_file = qv.ifEmpty(no_file('NO_QV')).first()
        RELEASE_REPORT(anno_ch, rel_gff, sym_summary, assembly, busco_summary, qv_file,
                       repeat_tbl, sets_report, mito_fasta.ifEmpty(no_file('NO_MITO_FASTA')).first(),
                       mito_check.ifEmpty(no_file('NO_MITO_CHECK')).first(), opt_file(params.release_notes, 'NO_NOTES'),
                       RELEASE_METHODS.out.methods)
        RELEASE_MQC_TABLE(anno_ch)

        def final_proteins = anno_ch.map { d -> file("${d}/annotate_results/*.proteins.fa") }.map { l -> l instanceof List ? l[0] : l }
        def stats_json     = anno_ch.map { d -> file("${d}/annotate_results/*.stats.json") }.map { l -> l instanceof List ? l[0] : l }
        ANNOTATION_QC_GATE(rel_gff, refined_gff, refine_report, protein_identity, sets_report, final_proteins,
            refined_proteins, sets.filter { f -> f.name == 'id_map.tsv' }.first(), gt_reports.collect(), stats_json,
            RELEASE_REPORT.out.summary, RELEASE_REPORT.out.readme, RELEASE_MQC_TABLE.out.table, SCAFFOLD_NAME_MAP.out.map,
            assembly, repeat_lib, classify_report, screen_files.collect().ifEmpty([no_file('NO_SCREEN_FILES')]),
            TABLE2ASN.out.files, sbt)

        def qc_files = ANNOTATION_QC_GATE.out.tsv.mix(refine_report, protein_identity, CURATE_PRODUCTS_FINAL.out.curation,
                                                       screen_files, sym_summary, qv_file, busco_summary)
            .flatten().filter { f -> !f.name.startsWith('NO_') }
        RELEASE_ASSEMBLE(ANNOTATION_QC_GATE.out.tsv, assembly, SCAFFOLD_NAME_MAP.out.map, anno_ch, sets.collect(),
            repeat_files.collect(), TABLE2ASN.out.files, mito_fasta.mix(mito_gff, mito_check).collect().ifEmpty([]),
            qc_files.unique { f -> f.name }.collect(), RELEASE_REPORT.out.readme.mix(RELEASE_REPORT.out.summary).collect())

    emit:
        anno     = anno_ch
        summary  = RELEASE_REPORT.out.summary
        ncbi     = TABLE2ASN.out.files
        qc_gate  = ANNOTATION_QC_GATE.out.tsv
        release  = RELEASE_ASSEMBLE.out.dir
        versions = versions_ch
}
