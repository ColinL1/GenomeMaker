// ============================================================
// 07b_annotation_qc.nf
// Annotation QC and ncRNA -- the parts of BRAKER4 that braker.pl does NOT do.
//
// SCOPE NOTE: braker.pl itself is deliberately untouched. BRAKER4's genemark/,
// augustus_training/ and augustus_predict/ rule sets reimplement what braker.pl
// already does internally, so porting those would mean rewriting braker.pl. What
// is genuinely additive is the QC, ncRNA and transcript-evidence work below.
//
// ALREADY COVERED, deliberately NOT ported:
//   best_by_compleasm -- braker.pl does this natively when --busco_lineage is
//     passed, which this pipeline already does. Verified in the v3.1.1 --help text
//     and by finding /opt/TSEBRA/bin/best_by_compleasm.py inside the image.
//   run_trnascan / run_barrnap -- already in stages 07 and 06.
//   run_masking -- stage 06 does this; repeats are a prerequisite here.
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: COMPLEASM
// Fast BUSCO-equivalent completeness on the genome.
//
// Faster than BUSCO on a 600 Mb genome and reports single/duplicated/fragmented
// in the same terms, so it cross-checks the BUSCO numbers rather than replacing
// them.
// ============================================================
process COMPLEASM {
    label 'annotation_qc_env'
    tag   "compleasm_${lineage}"

    conda 'bioconda::compleasm=0.2.9'
    container 'quay.io/biocontainers/compleasm:0.2.9--pyhdfd78af_0'
    publishDir "${params.outdir}/07_structural_annotation/qc/compleasm", mode: params.publish_mode

    input:
    path assembly
    val  lineage
    path compleasm_db

    output:
    path "compleasm_out/summary.txt", emit: summary
    path "compleasm_out",             emit: outdir
    path "versions.yml", emit: versions

    script:
    // compleasm keeps its lineage data in a download directory; pointing it at a
    // pre-populated one keeps this offline, exactly as BUSCO's --offline does.
    def db_arg = (compleasm_db && !compleasm_db.baseName.startsWith('NO_')) ? "-L ${compleasm_db}" : ""
    """
    compleasm run \\
      -a ${assembly} \\
      -o compleasm_out \\
      -l ${lineage} \\
      -t ${task.cpus} \\
      ${db_arg}

    if [ ! -s compleasm_out/summary.txt ]; then
        echo "ERROR: compleasm produced no summary -- check that lineage '${lineage}' is available" >&2
        echo "       in ${compleasm_db}, or allow it to download by unsetting params.compleasm_db." >&2
        exit 1
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        compleasm: "\$(compleasm --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p compleasm_out
    printf 'S:95.00%%, 950\\nD:1.00%%, 10\\nF:2.00%%, 20\\nM:2.00%%, 20\\nN:1000\\n' > compleasm_out/summary.txt
    echo "COMPLEASM stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        compleasm: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: OMARK
// Proteome completeness AND contamination, from OMAmer placements.
//
// WHY THIS MATTERS HERE MORE THAN BUSCO: braker.pl is run with --busco_lineage,
// which uses compleasm to select the TSEBRA mode that MAXIMISES BUSCO. The
// resulting BUSCO score is therefore optimised against, not independent of, the
// metric -- it cannot be reported as an unbiased completeness figure. OMArk is
// independent of that loop, so it is the honest headline number.
//
// OMArk also reports contamination, which is directly relevant to a coral
// holobiont where symbiont and bacterial sequence is expected.
// ============================================================
process OMARK {
    label 'annotation_qc_env'
    label 'res_anno_qc'
    tag   "omark"

    conda 'bioconda::omark=0.5.0'
    container 'quay.io/biocontainers/omark:0.5.0--pyhdfd78af_0'
    publishDir "${params.outdir}/07_structural_annotation/qc/omark", mode: params.publish_mode

    input:
    path proteins
    path gff3        // gene -> isoform structure, for OMArk's isoform file
    path omamer_db

    output:
    path "omark_out", emit: outdir
    path "omark_out/*.sum", emit: summary, optional: true
    path "versions.yml", emit: versions

    script:
    """
    # ete3 (used by OMArk for the NCBI taxonomy) creates ~/.etetoolkit and
    # downloads taxdump there; HOME is read-only in the container.
    export HOME=\$PWD

    # OMAmer placement first, then OMArk reads those placements.
    omamer search \\
      --db ${omamer_db} \\
      --query ${proteins} \\
      --out proteins.omamer \\
      --nthreads ${task.cpus}

    # OMArk's isoform ("splice") file: one line per gene, its transcript IDs
    # separated by ';'. Without it every isoform is scored as a separate protein
    # and alternative isoforms are counted as duplicated HOGs (36% "unexpected
    # duplicates" on this proteome vs 16% with it).
    awk -F'\\t' '!/^#/ && \$3=="mRNA" {
        match(\$9, /ID=[^;]+/);     id = substr(\$9, RSTART+3, RLENGTH-3)
        match(\$9, /Parent=[^;]+/); g  = substr(\$9, RSTART+7, RLENGTH-7)
        iso[g] = iso[g] (iso[g] ? ";" : "") id
    } END { for (g in iso) print iso[g] }' ${gff3} > proteins.splice

    mkdir -p omark_out
    omark \\
      -f proteins.omamer \\
      -d ${omamer_db} \\
      -i proteins.splice \\
      -o omark_out

    if ! ls omark_out/*.sum >/dev/null 2>&1; then
        echo "ERROR: OMArk produced no .sum report -- check the OMAmer database at ${omamer_db}" >&2
        exit 1
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omark: "\$(omark --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p omark_out
    printf 'Number of conserved HOGs: 1000\\nSingle:950\\nDuplicated:20\\nMissing:30\\n' > omark_out/stub.sum
    echo "OMARK stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        omark: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: INFERNAL_CMSCAN
// Rfam ncRNA annotation (rRNA, snRNA, snoRNA, miRNA, riboswitches...).
//
// This is a genuine gap rather than a refinement: the pipeline currently
// annotates tRNAs (tRNAscan-SE) and rRNAs (barrnap) and nothing else, so every
// other ncRNA class is simply absent from the final GFF3.
//
// -Z is REQUIRED for correct E-values: cmscan computes them against a stated
// search-space size in megabases, so it must be told the genome size (both
// strands) or the E-values are meaningless.
// ============================================================
process INFERNAL_CMSCAN {
    label 'annotation_qc_env'
    tag   "cmscan"

    conda 'bioconda::infernal=1.1.5'
    container 'quay.io/biocontainers/infernal:1.1.5--pl5321h031d066_3'
    publishDir "${params.outdir}/07_structural_annotation/ncrna", mode: params.publish_mode

    input:
    path assembly
    path rfam_dir      // WHOLE directory: cmscan needs Rfam.cm's cmpress indices
                       // (.i1f/.i1i/.i1m/.i1p) staged alongside the .cm itself
    val  rfam_cm_name
    val  rfam_clanin_name

    output:
    path "cmscan.tblout", emit: tblout
    path "cmscan.out",    emit: report
    path "versions.yml", emit: versions

    script:
    def clan_arg = rfam_clanin_name ? "--clanin ${rfam_dir}/${rfam_clanin_name}" : ""
    """
    # Search space in Mb, counting BOTH strands, as Infernal's manual requires.
    GSIZE=\$(grep -v '^>' ${assembly} | tr -d '\\n' | wc -c)
    Z=\$(awk -v g="\$GSIZE" 'BEGIN{printf "%.6f", g*2/1000000}')
    echo "genome 2x\$GSIZE bp -> -Z \$Z Mb"

    cmscan \\
      --cpu ${task.cpus} \\
      --rfam --cut_ga --nohmmonly \\
      --tblout cmscan.tblout \\
      --fmt 2 \\
      -Z \$Z \\
      ${clan_arg} \\
      ${rfam_dir}/${rfam_cm_name} ${assembly} > cmscan.out
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        infernal: "\$(cmscan -h 2>&1 | grep -oim1 'INFERNAL [0-9.]*' | sed 's/INFERNAL //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    printf '#idx target_name accession query_name\\n' > cmscan.tblout
    touch cmscan.out
    echo "INFERNAL_CMSCAN stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        infernal: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: STRINGTIE_ASSEMBLE
// Reference-guided transcript assembly from the merged RNA-seq BAM.
//
// SCOPE, stated honestly: this produces UTR-bearing transcript models as
// EVIDENCE and reports transcript-level statistics. It does NOT rewrite the
// BRAKER gene models. Faithfully transferring UTRs onto predicted CDS requires
// GUSHR (which is what BRAKER4's add_utr rule uses) and GUSHR is not present in
// the braker3 image or packaged in bioconda. Doing the transfer approximately
// would corrupt gene models, which is worse than not having UTRs, so the
// transfer is left as a separate deliberate step.
// ============================================================
process STRINGTIE_ASSEMBLE {
    label 'annotation_qc_env'
    label 'res_anno_qc'
    tag   "stringtie"

    conda 'bioconda::stringtie=3.0.3'
    container 'quay.io/biocontainers/stringtie:3.0.3--h29c0135_0'
    publishDir "${params.outdir}/07_structural_annotation/stringtie", mode: params.publish_mode

    input:
    path bam
    path guide_gff

    output:
    path "stringtie.gtf",       emit: gtf
    path "stringtie_stats.txt", emit: stats
    path "versions.yml", emit: versions

    script:
    // -G guides assembly with the existing models so transcript IDs stay
    // relatable to the gene set; --conservative avoids inventing structures that
    // the RNA-seq does not clearly support.
    """
    stringtie ${bam} \\
      -G ${guide_gff} \\
      -p ${task.cpus} \\
      --conservative \\
      -o stringtie.gtf

    {
      echo "# StringTie reference-guided transcript assembly"
      printf 'transcripts\\t%s\\n' "\$(awk '\$3=="transcript"' stringtie.gtf | wc -l)"
      printf 'exons\\t%s\\n'       "\$(awk '\$3=="exon"' stringtie.gtf | wc -l)"
      printf 'with_utr\\t%s\\n'    "\$(grep -c 'UTR' stringtie.gtf || true)"
    } > stringtie_stats.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: "\$(stringtie --version 2>&1 | head -1 || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    printf '# stub\\n' > stringtie.gtf
    printf 'transcripts\\t0\\n' > stringtie_stats.txt
    echo "STRINGTIE_ASSEMBLE stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stringtie: "stub"
    END_VERSIONS
    """
}
