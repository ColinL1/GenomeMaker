// ============================================================
// subworkflows/03_assembly_filtering.nf
// BlobToolKit-based assembly filtering -- WORKFLOW ONLY
// Process definitions live in modules/03_assembly_filtering.nf
//
// Runs on the unpurged contigs (after the mito-contig removal and polishing), so the host and
// symbiont bins are decided by taxonomy before purge_dups' coverage cutoffs see them:
//   ASSEMBLY_FILTERING(contigs, porechopped_ch, busco_table,
//                      blast_db_path, blast_taxdump_path, diamond_db_path)
// with two named outputs downstream code depends on:
//   .out.filtered_assembly  -- host-only fasta, feeds PURGE_DUPS
//   .out.symbiont_assembly  -- symbiont-only fasta
//
// Categories/taxon names to keep are NOT passed in by main.nf, so they're
// built here from params (with defaults) rather than as a take() input.
// ASSUMPTION: add these to nextflow.config if you want non-default taxa:
//   params.host_taxon_name     = 'Cnidaria'
//   params.symbiont_taxon_name = 'Symbiodiniaceae'
// If run_blobtools=false, symbiont_assembly emits an empty channel since
// PASSTHROUGH never separates a symbiont fraction.
// ============================================================

nextflow.enable.dsl = 2

include { WRITE_META_YAML       } from '../modules/03_assembly_filtering.nf'
include { PREPARE_BLAST_QUERY   } from '../modules/03_assembly_filtering.nf'
include { BLASTN_NT             } from '../modules/03_assembly_filtering.nf'
include { DIAMOND_BLASTX        } from '../modules/03_assembly_filtering.nf'
include { MINIMAP2_SORT         } from '../modules/03_assembly_filtering.nf'
include { BLOBTOOLS_CREATE      } from '../modules/03_assembly_filtering.nf'
include { BLOBTOOLS_ADD_HITS    } from '../modules/03_assembly_filtering.nf'
include { BLOBTOOLS_ADD_COV     } from '../modules/03_assembly_filtering.nf'
include { BLOBTOOLS_ADD_BUSCO   } from '../modules/03_assembly_filtering.nf'
include { BLOBTOOLS_APPLY_FILTER } from '../modules/03_assembly_filtering.nf'
include { BLOBTOOLS_STATIC_PLOTS } from '../modules/03_assembly_filtering.nf'
include { PASSTHROUGH           } from '../modules/03_assembly_filtering.nf'

workflow ASSEMBLY_FILTERING {

    take:
        assembly
        reads
        busco_table
        blast_db_dir   // directory bind-mounted into the container so BLASTN_NT/BLOBTOOLS_CREATE can see it
        taxdump_dir    // NCBI new_taxdump/ sibling directory to blast_db_dir
        diamond_db     // UniProt reference proteomes .dmnd file

    main:
        def versions_ch = channel.empty()
        def report_ch   = channel.empty()   // files for the MultiQC report

        if (!params.run_blobtools.toString().toBoolean()) {
            // Skip filtering, pass assembly through
            PASSTHROUGH(assembly)
            filtered  = PASSTHROUGH.out.filtered_assembly
            symbiont  = channel.empty()
            cov_stats = channel.empty()
            versions_ch = versions_ch.mix(PASSTHROUGH.out.versions)
        } else {
            // -- evidence generation, one process per tool --
            // BLASTN and DIAMOND are deliberately fed DIFFERENT query sets,
            // because the two tools have opposite cost models:
            //
            //   BLASTN  scales with query x database, so splitting the work
            //           across tasks is a real win -- but only if each chunk
            //           holds many records, otherwise -mt_mode 1 cannot
            //           parallelise. PREPARE_BLAST_QUERY windows the handful
            //           of chromosome-scale scaffolds so every chunk is a
            //           uniform bundle of many queries.
            //
            //   DIAMOND is dominated by building the seed index over the
            //           32.0 Gletter reference, which is paid in full on
            //           every invocation no matter how small the query is.
            //           Chunking therefore multiplies the cost by the chunk
            //           count for no benefit, so it gets the whole assembly
            //           in a single call. See modules/03_assembly_filtering.nf
            //           for the measurements behind both decisions.
            PREPARE_BLAST_QUERY(assembly)

            // TWO independent chunkings of the SAME windowed query, because the two
            // tools have different optima -- both MEASURED, see the cost models in
            // modules/03_assembly_filtering.nf:
            //   blastn  : 5 MB  -- many small records so -mt_mode 1 can split the
            //             query across threads; startup is only ~5 min.
            //   diamond : 40 MB -- startup is ~12.6 min and the search itself is
            //             nearly free (4 s/MB), so bigger chunks amortise far
            //             better. 5 MB would be 4x slower per MB.
            def blast_query = PREPARE_BLAST_QUERY.out.query

            blastn_chunks  = blast_query.splitFasta(size: params.blastn_chunk_size,  file: true)
            diamond_chunks = blast_query.splitFasta(size: params.diamond_chunk_size, file: true)

            BLASTN_NT(blastn_chunks, blast_db_dir)
            DIAMOND_BLASTX(diamond_chunks, diamond_db)

            blastn_merged = BLASTN_NT.out.blastn_out.collectFile(name: 'blastn.out')
            blastx_merged = DIAMOND_BLASTX.out.blastx_out.collectFile(name: 'blastx.out')

            // One process: minimap2 piped straight into samtools sort, so no
            // uncompressed SAM is ever written to the work dir.
            MINIMAP2_SORT(assembly, reads)

            // -- BlobDir build, chained one `blobtools add` per step --
            // meta.yaml is generated by a task, not written into workDir at parse time.
            WRITE_META_YAML()

            BLOBTOOLS_CREATE(assembly, WRITE_META_YAML.out.meta_yaml, params.taxon_id, taxdump_dir)

            BLOBTOOLS_ADD_HITS(
                BLOBTOOLS_CREATE.out.blobdir,
                blastn_merged,
                blastx_merged,
                taxdump_dir
            )

            BLOBTOOLS_ADD_COV(
                BLOBTOOLS_ADD_HITS.out.blobdir,
                MINIMAP2_SORT.out.cov_bam
            )

            BLOBTOOLS_ADD_BUSCO(
                BLOBTOOLS_ADD_COV.out.blobdir,
                busco_table
            )

            blobdir_final = BLOBTOOLS_ADD_BUSCO.out.blobdir

            // -- filtering modes (the BlobDir above is always produced/published) --
            //   1. params.filter_json set : HOST = manual selection (one JSON list
            //      exported from `blobtools view --remote btk`, BlobtoolKit Viewer)
            //   2. params.auto_filter     : automatic best-hit filter (default on)
            //   3. neither                : no filtering, the assembly goes on as-is
            // Auto filter: HOST keeps Cnidaria (+ no-hit), SYMBIONT keeps
            // Symbiodiniaceae. Each category carries its OWN field because the taxa
            // sit at different ranks (phylum vs FAMILY). A manual JSON replaces the
            // host selection only; the symbiont set stays automatic if auto_filter
            // is on. No intermediate JSON is generated. `[]` stages no file.
            // A CLI `--auto_filter false` arrives as the STRING "false", which is
            // truthy in Groovy, so parse it explicitly.
            def auto_filter = params.auto_filter.toString().toBoolean()
            def manual_json = params.filter_json ? file(params.filter_json, checkIfExists: true) : []
            def cats = []
            if (params.filter_json || auto_filter) {
                cats << ['host',     params.host_taxon_rank     ?: 'bestsumorder_phylum',
                                     params.host_taxon_name     ?: 'Cnidaria',
                                     manual_json]
            }
            if (auto_filter) {
                cats << ['symbiont', params.symbiont_taxon_rank ?: 'bestsumorder_family',
                                     params.symbiont_taxon_name ?: 'Symbiodiniaceae',
                                     []]
            }

            if (cats) {
                apply_input = channel.fromList(cats)
                    .combine(assembly)
                    .combine(blobdir_final)

                BLOBTOOLS_APPLY_FILTER(apply_input)

                filtered = BLOBTOOLS_APPLY_FILTER.out.filtered_assembly
                    .filter { category, path -> category == 'host' }
                    .map    { category, path -> path }

                symbiont = BLOBTOOLS_APPLY_FILTER.out.filtered_assembly
                    .filter { category, path -> category == 'symbiont' }
                    .map    { category, path -> path }

                versions_ch = versions_ch.mix(BLOBTOOLS_APPLY_FILTER.out.versions.first())
                report_ch = report_ch.mix(
                    BLOBTOOLS_APPLY_FILTER.out.summary.map { category, p -> p },
                    BLOBTOOLS_APPLY_FILTER.out.scaffolds)
            } else {
                // mode 3: BlobDir only, assembly continues unfiltered
                filtered = assembly
                symbiont = channel.empty()
            }

            // Static plots of the FINAL OUTPUT of this stage, once (params.blobtools_plots,
            // default on): the filtered host BlobDir when a filter ran (automatic or
            // manual JSON), otherwise the unfiltered final BlobDir.
            // Parsed explicitly: a CLI `--blobtools_plots false` arrives as a truthy string.
            if (params.blobtools_plots.toString().toBoolean()) {
                plot_in = cats
                    ? BLOBTOOLS_APPLY_FILTER.out.filtered_blobdir
                        .filter { category, dir -> category == 'host' }
                        .map    { category, dir -> dir }
                    : blobdir_final
                BLOBTOOLS_STATIC_PLOTS(plot_in)
                versions_ch = versions_ch.mix(BLOBTOOLS_STATIC_PLOTS.out.versions)
                report_ch   = report_ch.mix(BLOBTOOLS_STATIC_PLOTS.out.plots.flatten())
            }

            cov_stats = MINIMAP2_SORT.out.flagstat

            versions_ch = versions_ch.mix(
                WRITE_META_YAML.out.versions, PREPARE_BLAST_QUERY.out.versions,
                BLASTN_NT.out.versions.first(), DIAMOND_BLASTX.out.versions.first(),
                MINIMAP2_SORT.out.versions, BLOBTOOLS_CREATE.out.versions,
                BLOBTOOLS_ADD_HITS.out.versions, BLOBTOOLS_ADD_COV.out.versions,
                BLOBTOOLS_ADD_BUSCO.out.versions)
        }

    emit:
        filtered_assembly = filtered
        symbiont_assembly = symbiont
        coverage_flagstat = cov_stats
        report_files      = report_ch       // blobtools filter summaries + static plots
        versions          = versions_ch
}
