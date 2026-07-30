// ============================================================
// 08_functional_annotation.nf
// Functional annotation pipeline:
// Phobius → InterProScan (funannotate iprscan) → eggNOG-mapper → funannotate annotate
// ============================================================

nextflow.enable.dsl = 2

include { PHOBIUS_RUN ; IPRSCAN_RUN ; EGGNOG_MAPPER ; FUNANNOTATE_ANNOTATE ; FUNANNOTATE_STATS } from '../modules/08_functional_annotation.nf'

// ============================================================
// Workflow: FUNCTIONAL_ANNOTATION
// Phobius → InterProScan → eggNOG → funannotate annotate → stats
// ============================================================
workflow FUNCTIONAL_ANNOTATION {

    take:
        assembly
        gff3
        proteins
        funannotate_db

    main:
        PHOBIUS_RUN(proteins)
        IPRSCAN_RUN(proteins)
        EGGNOG_MAPPER(proteins, funannotate_db)
        FUNANNOTATE_ANNOTATE(
            gff3,
            assembly,
            EGGNOG_MAPPER.out.eggnog_anno,
            IPRSCAN_RUN.out.iprscan_xml,
            PHOBIUS_RUN.out.phobius_results
        )
        FUNANNOTATE_STATS(assembly, gff3)

    emit:
        anno_dir   = FUNANNOTATE_ANNOTATE.out.anno_dir
        stats_json = FUNANNOTATE_STATS.out.stats_json
}
