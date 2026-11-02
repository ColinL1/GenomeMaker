// ============================================================
// subworkflows/03_assembly_filtering.nf
// BlobToolKit host/symbiont split on the unpurged contigs; processes in
// modules/03_assembly_filtering.nf. Without run_blobtools there is no symbiont bin.
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
        blast_db_dir
        taxdump_dir   
        diamond_db

    main:
        def versions_ch = channel.empty()
        def report_ch   = channel.empty()

        if (!params.run_blobtools.toString().toBoolean()) {
            PASSTHROUGH(assembly)
            filtered  = PASSTHROUGH.out.filtered_assembly
            symbiont  = channel.empty()
            cov_stats = channel.empty()
            versions_ch = versions_ch.mix(PASSTHROUGH.out.versions)
        } else {
            PREPARE_BLAST_QUERY(assembly)

            // separate chunk sizes: diamond's startup dominates, so its chunks are larger
            def blast_query = PREPARE_BLAST_QUERY.out.query

            blastn_chunks  = blast_query.splitFasta(size: params.blastn_chunk_size,  file: true)
            diamond_chunks = blast_query.splitFasta(size: params.diamond_chunk_size, file: true)

            BLASTN_NT(blastn_chunks, blast_db_dir)
            DIAMOND_BLASTX(diamond_chunks, diamond_db)

            blastn_merged = BLASTN_NT.out.blastn_out.collectFile(name: 'blastn.out')
            blastx_merged = DIAMOND_BLASTX.out.blastx_out.collectFile(name: 'blastx.out')

            MINIMAP2_SORT(assembly, reads)

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

            // Filter modes: filter_json (manual host selection), auto_filter (best hit), or
            // neither (BlobDir only). CLI booleans arrive as strings, hence toBoolean().
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

                // first(): the host bin as a value channel, like the other branches
                filtered = BLOBTOOLS_APPLY_FILTER.out.filtered_assembly
                    .filter { category, path -> category == 'host' }
                    .map    { category, path -> path }
                    .first()

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

            // Static plots of the filtered host BlobDir, else the unfiltered one
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
        report_files      = report_ch       
        versions          = versions_ch
}
