// ============================================================
// 07b_annotation_qc.nf
// Annotation QC and ncRNA: compleasm, OMArk, Infernal/Rfam, StringTie.
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: COMPLEASM
// Fast BUSCO-equivalent completeness on the genome.
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
    // pre-populated lineage dir keeps compleasm offline
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
// OMArk is independent of BRAKER's BUSCO-driven TSEBRA mode choice, and also
// reports contamination.
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
    # ete3 writes ~/.etetoolkit; HOME is read-only in the container
    export HOME=\$PWD

    omamer search \\
      --db ${omamer_db} \\
      --query ${proteins} \\
      --out proteins.omamer \\
      --nthreads ${task.cpus}

    # OMArk's isoform ("splice") file: one line per gene, transcript IDs separated
    # by ';'. Without it isoforms count as duplicated HOGs.
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
// -Z: search space in Mb (both strands), needed for correct E-values.
// ============================================================
process INFERNAL_CMSCAN {
    label 'annotation_qc_env'
    tag   "cmscan"

    conda 'bioconda::infernal=1.1.5'
    container 'quay.io/biocontainers/infernal:1.1.5--pl5321h031d066_3'
    publishDir "${params.outdir}/07_structural_annotation/ncrna", mode: params.publish_mode

    input:
    path assembly
    path rfam_dir      // whole directory: cmscan needs the cmpress indices next to Rfam.cm
    val  rfam_cm_name
    val  rfam_clanin_name

    output:
    path "cmscan.tblout", emit: tblout
    path "cmscan.out",    emit: report
    path "versions.yml", emit: versions

    script:
    def clan_arg = rfam_clanin_name ? "--clanin ${rfam_dir}/${rfam_clanin_name}" : ""
    """
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
// Transcript models are evidence and statistics only; the gene models are not changed.
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
