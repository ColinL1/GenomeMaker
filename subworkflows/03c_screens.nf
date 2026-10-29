// ============================================================
// 03c_screens.nf (subworkflow)
// NCBI contamination screens on the host-filtered contigs, before scaffolding:
// FCS-adaptor + VecScreen + FCS-GX, review of the unclear hits, then the actions
// are applied. The sequence is not changed after this step.
// ============================================================

nextflow.enable.dsl = 2

include { FCS_ADAPTOR ; VECSCREEN ; VECSCREEN_CLASSIFY ; FCS_GX ; SCREEN_REVIEW ; APPLY_SCREEN_ACTIONS } from '../modules/03c_screens.nf'

workflow SCREEN_ASSEMBLY {

    take:
        assembly        // value channel: host-filtered assembly

    main:
        def in_ch = assembly.map { f -> ['contigs', f] }
        def decisions = params.screen_decisions
            ? file(params.screen_decisions, checkIfExists: true)
            : file("${projectDir}/assets/NO_SCREEN_DECISIONS", checkIfExists: true)

        FCS_ADAPTOR(in_ch, file(params.fcs_adaptor_image, checkIfExists: true))
        VECSCREEN(in_ch, file(params.univec, checkIfExists: true))
        VECSCREEN_CLASSIFY(VECSCREEN.out.raw)
        FCS_GX(in_ch, file(params.fcs_gx_db, checkIfExists: true), file(params.fcs_gx_image, checkIfExists: true))

        SCREEN_REVIEW(FCS_GX.out.report.join(VECSCREEN_CLASSIFY.out.regions), decisions)

        // [sname, fasta, [gx report, adaptor report], review_excluded, review_mask, screen_review]
        def apply_in = in_ch
            .join(FCS_GX.out.report.join(FCS_ADAPTOR.out.report).map { s, gx, ad -> [s, [gx, ad]] })
            .join(SCREEN_REVIEW.out.actions)
        APPLY_SCREEN_ACTIONS(apply_in)

    emit:
        screened = APPLY_SCREEN_ACTIONS.out.fasta
        qc_files = APPLY_SCREEN_ACTIONS.out.log.mix(SCREEN_REVIEW.out.review)   // screen_actions.tsv, screen_review.tsv
        versions = FCS_ADAPTOR.out.versions.mix(VECSCREEN.out.versions, FCS_GX.out.versions)
}
