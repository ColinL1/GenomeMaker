// ============================================================
// 08_functional_annotation.nf
// eggNOG-mapper → (InterProScan) → (Phobius) → funannotate annotate → stats
//
// InterProScan runs when its data directory is configured, Phobius with params.run_phobius.
// ============================================================

nextflow.enable.dsl = 2

include { PHOBIUS_RUN ; INTERPROSCAN ; EGGNOG_MAPPER ; FUNANNOTATE_ANNOTATE } from '../modules/08_functional_annotation.nf'

workflow FUNCTIONAL_ANNOTATION {

    take:
        assembly
        gff3
        proteins
        funannotate_db
        eggnog_db

    main:
        // distinct placeholder names: one shared NO_FILE causes an input file name collision
        def no_ipr = channel.value(file("${projectDir}/assets/NO_IPRSCAN.xml", checkIfExists: true))
        def no_pho = channel.value(file("${projectDir}/assets/NO_PHOBIUS.txt", checkIfExists: true))

        def versions_ch = channel.empty()

        EGGNOG_MAPPER(proteins, eggnog_db)

        // --- optional: InterProScan ---
        def ipr_ch = no_ipr
        // local install (interproscan_dir) or container + matching data (interproscan_data)
        def ipr_local     = params.interproscan_dir && file("${params.interproscan_dir}/interproscan.sh").exists()
        def ipr_container = !params.interproscan_dir && params.interproscan_data && file(params.interproscan_data).exists()
        if (ipr_local || ipr_container) {
            INTERPROSCAN(proteins, channel.value(file(ipr_local ? params.interproscan_dir : params.interproscan_data)))
            ipr_ch = INTERPROSCAN.out.iprscan_xml
            versions_ch = versions_ch.mix(INTERPROSCAN.out.versions)
        } else {
            log.warn "InterProScan SKIPPED (interproscan_dir=${params.interproscan_dir}) -- " +
                        "funannotate will run without --iprscan, so the annotation will carry NO " +
                        "InterPro domain assignments and none of the GO terms derived from them"
        }

        // --- optional: Phobius ---
        def pho_ch = no_pho
        if (params.run_phobius.toString().toBoolean()) {
            PHOBIUS_RUN(proteins, channel.value(file(params.phobius_dir, checkIfExists: true)))
            pho_ch = PHOBIUS_RUN.out.phobius_results
            versions_ch = versions_ch.mix(PHOBIUS_RUN.out.versions)
        } else {
            log.warn "Phobius SKIPPED (run_phobius=false) -- no signal-peptide or " +
                        "transmembrane predictions in the annotation"
        }

        FUNANNOTATE_ANNOTATE(
            gff3,
            assembly,
            EGGNOG_MAPPER.out.eggnog_anno,
            ipr_ch,
            pho_ch,
            funannotate_db
        )
        // statistics from funannotate's own stats.json
        def stats_json_ch = FUNANNOTATE_ANNOTATE.out.anno_dir
            .map { d -> file("${d}/annotate_results/*.stats.json") }
            .map { l -> l instanceof List ? l[0] : l }

        versions_ch = versions_ch.mix(
            EGGNOG_MAPPER.out.versions, FUNANNOTATE_ANNOTATE.out.versions)

    emit:
        anno_dir   = FUNANNOTATE_ANNOTATE.out.anno_dir
        fixedproducts = FUNANNOTATE_ANNOTATE.out.fixedproducts
        stats_json = stats_json_ch
        eggnog     = EGGNOG_MAPPER.out.eggnog_anno
        iprscan    = ipr_ch
        phobius    = pho_ch
        versions   = versions_ch
}
