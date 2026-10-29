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
include { linkage_evidence }               from '../modules/10_release.nf'

// Methods template of the run: one line per step that ran, with the settings it used.
// {{tool}} becomes the version the run recorded and {{=result}} a sentence read from the
// decision files (bin/methods.py, process METHODS); no result is written here.
// run: [input_assembly, mito, structural, functional, release] -- what main.nf decided.
def base_name(f) {
    return file(f.toString()).name
}

def methods_template(Map run) {
    def lines = []
    def porechop = params.run_porechop.toString().toBoolean() ? 'adapter-trimmed with Porechop {{porechop}}, ' : ''
    def contam = (params.symbiont_refs && params.run_symbiont_filter.toString().toBoolean())
        ? '; reads matching the symbiont reference removed (chopper --contam)' : ''
    lines << "- Reads: ONT reads ${porechop}split with chopper {{chopper}} into an assembly set (q >= 3, >= 1 kb) " +
             "and a polishing set (q >= 5, >= 500 bp)${contam}."
    if (run.input_assembly) {
        lines << "- Assembly: supplied assembly `${base_name(params.input_assembly)}`."
    } else {
        def names = [necat: 'NECAT {{necat}}', flye: 'Flye {{flye}}', hifiasm: 'hifiasm {{hifiasm}}', raven: 'Raven {{raven}}']
        def asm = params.assemblers.collect { a -> names[a.toString()] ?: a.toString() }.join(', ')
        def size = params.genome_size
            ? "${params.genome_size} bp (given)"
            : "GenomeScope2 {{genomescope2}} estimate (meryl {{meryl}}, k = ${params.kmer_size})"
        lines << "- Assembly: ${asm}; expected genome size ${size}. Assemblies ranked by BUSCO {{busco}} " +
                 "(${params.busco_lineage}) completeness, size outliers (more than ${Math.round((params.assembly_size_tolerance as double) * 100)}% from the expected size) last; " +
                 "within ${params.busco_tie_band} pp of the best completeness the assembly closest to the expected size wins: {{=assembly}}."
    }
    lines << (run.mito
        ? "- Mitochondrial contigs: MitoHiFi {{mitohifi}} (contigs mode, reference `${base_name(params.mito_ref)}`) on the selected assembly; " +
          "contigs it lists, up to ${params.mito_contig_max_ratio} x the reference length, removed from the nuclear assembly: {{=mito}}."
        : "- Mitochondrial contigs: not checked (run_mitogenome=false).")
    lines << (params.run_polish.toString().toBoolean()
        ? "- Polishing: Medaka {{medaka}} with the polishing reads, on the contigs; {{=polish}}."
        : "- Polishing: not run (run_polish=false).")
    if (!params.run_blobtools.toString().toBoolean()) {
        lines << "- Host/symbiont split: not run (run_blobtools=false)."
    } else {
        def sel = params.filter_json
            ? "host contigs from the manual selection `${base_name(params.filter_json)}`"
            : "host = ${params.host_taxon_rank} ${params.host_taxon_name} and no-hit contigs"
        def sym = params.auto_filter.toString().toBoolean()
            ? "; symbiont bin = ${params.symbiont_taxon_rank} ${params.symbiont_taxon_name}" : ''
        lines << "- Host/symbiont split: BlobToolKit {{blobtoolkit}} (BLASTn {{blast}} against `${base_name(params.blast_db)}`, DIAMOND {{diamond}}, " +
                 "ONT coverage, BUSCO) on the unpurged contigs; ${sel}${sym}; other contigs removed."
    }
    if (params.run_purge_dups.toString().toBoolean()) {
        def cut = params.purge_dups_transition != null
            ? "calcuts transition ${params.purge_dups_transition} (lower ${params.purge_dups_lower ?: 'automatic'}, upper ${params.purge_dups_upper ?: 'automatic'})"
            : 'automatic calcuts cutoffs'
        lines << "- Haplotig purging: purge_dups {{purge_dups}} on the host contigs, ${cut}, `get_seqs -e`; {{=purge}}."
    } else {
        lines << "- Haplotig purging: not run (run_purge_dups=false)."
    }
    lines << "- Length filter: sequences shorter than ${params.min_contig_length} bp removed: {{=length}}."
    if (params.run_screens.toString().toBoolean()) {
        lines << "- Contamination screens: FCS-adaptor {{fcs-adaptor}} (--euk), VecScreen (UniVec_Core, BLAST {{blast}}) and FCS-GX {{fcs-gx}} " +
                 "(tax-id ${params.taxon_id}) on the contigs; unclear hits reviewed (`screen_review.tsv`); {{=screen}}. " +
                 "Later steps add scaffold gaps and rename sequences only."
    } else {
        lines << "- Contamination screens: not run (run_screens=false)."
    }
    if (params.scaffold_reference) {
        lines << "- Scaffolding: RagTag {{ragtag}} `scaffold` against `${base_name(params.scaffold_reference)}`" +
                 (params.ragtag_infer_gaps.toString().toBoolean() ? ' (gap sizes inferred)' : ' (100 bp gaps of unknown size)') + '.'
    }
    lines << "- Naming: sequences sorted by length and renamed ${params.scaffold_basename}_1..N (funannotate sort {{funannotate}})."
    if (run.mito) {
        def route = params.illumina_dna_r1
            ? 'GetOrganelle {{getorganelle}} on the Illumina DNA reads'
            : 'the MitoHiFi contig polished with Medaka {{medaka}} using the mito-mapped ONT reads'
        lines << "- Mitogenome: ${route}; MITOS2 {{mitos}} annotation" +
                 (params.mito_annot_ref_fasta && !(params.mito_annot_ref_fasta.toString() in ['false', 'null', ''])
                     ? ", Liftoff {{liftoff}} from `${base_name(params.mito_annot_ref_fasta)}`" : '') + '.'
    }
    if (params.run_repeats.toString().toBoolean()) {
        lines << "- Repeats: RepeatModeler {{repeatmodeler}} + EDTA {{edta}} library; Unknown families classified with RepeatClassifier; " +
                 "RepeatMasker {{repeatmasker}} soft-masking with the classified library."
    }
    if (run.structural) {
        lines << "- Gene models: BRAKER4 {{braker4}} (RNA-seq + proteins), tRNAscan-SE {{trnascan-se}} high-confidence set, " +
                 "Rfam/Infernal {{infernal}} ncRNAs; UTRs from BRAKER's stringtie2utr, gene membership by a frame-aware " +
                 "CDS-overlap rule (>= ${params.refine_frame_min}), Rfam filters, evidence and haplotig tags" +
                 (params.related_proteomes ? ", models fused across assembly gaps split where >= ${params.gap_split_min_species} related proteomes support it" : '') +
                 (params.related_proteomes && params.isoform_filter.toString().toBoolean()
                     ? "; alternative isoforms kept only with evidence of their own (an intron with >= ${params.isoform_min_reads} RNA-seq reads, " +
                       "RNA-seq depth, or better cover of a related protein), one isoform kept per gene" : '') + '.'
    }
    if (run.functional) {
        lines << "- Functional annotation: eggNOG-mapper {{eggnog-mapper}}" +
                 (params.interproscan_data || params.interproscan_dir ? ', InterProScan {{interproscan}}' : '') +
                 (params.run_phobius.toString().toBoolean() ? ', Phobius {{phobius}}' : '') + ", funannotate {{funannotate}}."
    }
    if (run.release) {
        lines << "- Release curation: product names checked against UniProt reviewed entries by gene symbol (UniProt release {{uniprot}}; " +
                 "`product_curation.tsv`); funannotate NCBI pass; table2asn {{table2asn}} with the organism '${params.ncbi_species}'" +
                 (params.locus_tag ? ", locus-tag prefix ${params.locus_tag}" : '') + "; gaps U / 100 bp / ${linkage_evidence()}; " +
                 "EC numbers and product names table2asn changed in a preflight run applied to the submitted annotation."
    }
    return lines.join('\n')
}

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
        custom_content    = MULTIQC_CUSTOM_CONTENT.out.mqc   // *_mqc tables
        multiqc_data      = MULTIQC_FINAL.out.data
        software_versions = COLLATE_VERSIONS.out.table
        pipeline_report   = GENERATE_SUMMARY.out.pipeline_report
}
