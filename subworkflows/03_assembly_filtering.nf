// ============================================================
// 03_assembly_filtering.nf
// BlobTools filtering or passthrough
//
// Params expected (add to nextflow.config):
//   params.skip_blobtools        = false
//   params.taxon_id              = <host NCBI taxid>       // anchors BlobDir lineage, NOT a filter
//   params.host_taxon_name       = 'Cnidaria'               // required keep-filter target
//   params.symbiont_taxon_name   = 'Dinophyceae'             // optional; set '' / null to disable
//
// Emits:
//   filtered_assembly  -> single path, HOST assembly only (same contract as before;
//                         main.nf's ASSEMBLY_POLISHING call needs no changes)
//   symbiont_assembly  -> single path, symbiont draft, or empty channel if disabled/skipped
//   filter_summary     -> per-category summary.json files, for the summary report
// ============================================================

nextflow.enable.dsl = 2

include { BLAST_AND_COV; GENERATE_TAXON_FILTER; BLOBTOOLS_APPLY_FILTER; PASSTHROUGH } from '../modules/03_assembly_filtering.nf'

workflow ASSEMBLY_FILTERING {

    take:
        assembly      // value channel, single path (chosen "best" assembly)
        reads         // ONT reads for coverage (porechopped_ch)
        busco_table   // single BUSCO full_table.tsv, anthozoa_odb12
        meta_yaml     // sample metadata YAML

    main:
        if (params.skip_blobtools) {
            PASSTHROUGH(assembly)
            filtered_assembly = PASSTHROUGH.out.filtered_assembly
            symbiont_assembly = channel.empty()
            filter_summary    = channel.empty()

        } else {
            BLAST_AND_COV(assembly, reads, busco_table, meta_yaml, params.taxon_id)

            // Build the (category, taxon_name) fan-out. Host is always run;
            // symbiont only runs if a target taxon is configured.
            def taxon_targets_ch = channel.of(['host', params.host_taxon_name])
            if (params.symbiont_taxon_name) {
                taxon_targets_ch = taxon_targets_ch.mix(
                    channel.of(['symbiont', params.symbiont_taxon_name])
                )
            }

            GENERATE_TAXON_FILTER(
                taxon_targets_ch.combine(BLAST_AND_COV.out.blobdir)
            )

            BLOBTOOLS_APPLY_FILTER(
                GENERATE_TAXON_FILTER.out.filter_json
                    .combine(assembly)
                    .combine(BLAST_AND_COV.out.blobdir)
            )

            // Split the (category, path) stream back into named outputs so
            // downstream consumers (e.g. ASSEMBLY_POLISHING in main.nf) don't
            // need to know about the category-tuple shape at all.
            filtered_assembly = BLOBTOOLS_APPLY_FILTER.out.filtered_assembly
                .filter { category, path -> category == 'host' }
                .map    { category, path -> path }

            symbiont_assembly = params.symbiont_taxon_name
                ? BLOBTOOLS_APPLY_FILTER.out.filtered_assembly
                    .filter { category, path -> category == 'symbiont' }
                    .map    { category, path -> path }
                : channel.empty()

            filter_summary = GENERATE_TAXON_FILTER.out.summary
                .map { category, path -> path }
        }

    emit:
        filtered_assembly
        symbiont_assembly
        filter_summary
}
