// ============================================================
// 04_assembly_polish.nf
// Assembly polishing: Racon (ONT) → Medaka (ONT) → funannotate clean/sort
// Final BUSCO on polished assembly
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// NOTE: a RACON_POLISH process lived here and has been REMOVED.
//
// It was also broken: it called minimap2 while its container was
// racon:1.5.0, which does not ship minimap2.
//
// Why it is gone rather than fixed:
//   * Medaka's models are trained on draft assemblies from particular
//     pipelines; feeding it a Racon-modified draft moves the input off that
//     training distribution. ONT's own medaka_consensus workflow dropped the
//     Racon pre-step for this reason.
//   * On R10.4.1 + Dorado SUP the reads are already ~Q20+, so there is little
//     for an overlap-consensus step to fix, and it can introduce errors in
//     repetitive regions.
//   * It cost a full extra minimap2 pass over 14.5 GB of reads.
//
// Merqury QV (MERQURY_QV below) is the measurement that would settle this
// empirically if it is ever revisited: git history has the original process.
// ============================================================

// ============================================================
// Process: MEDAKA_POLISH
// Medaka consensus polishing with ONT reads.
//
// Uses the `medaka_consensus` wrapper, which does its own mini_align ->
// sorted BAM internally. That is deliberate: the previous version hand-rolled
// minimap2 -> SAM -> samtools sort, writing a full uncompressed SAM to disk for
// a 600 Mb assembly against 14.5 GB of reads.
//
// The previous version also called `medaka.consensus`, which is not a command in
// any medaka release. Verified against medaka 2.2.2: the subcommands are
// {compress_bam, features, train, inference, smolecule, tandem,
// consensus_from_features, fastrle, sequence, vcf, tools} -- there is no
// `consensus` subcommand either (2.x renamed it to `inference` + `sequence`),
// so `medaka_consensus` is the correct entry point.
//
// MODEL SELECTION: medaka's model must match the basecaller, and a mismatch
// silently degrades consensus accuracy -- so this never guesses. The old
// fallback, `medaka list_models | head -1 | awk '{print $1}'`, picked the
// header line of that listing rather than a model. Order of preference:
//   1. params.medaka_model, if set
//   2. the basecaller model from the ONT read headers (basecall_model_version_id=
//      or the Dorado RG:Z: tag), validated via
//      `medaka tools resolve_model`
//   3. hard failure with an actionable message
// check for Dosnt't medaka now have a "medaka tools resolve_model --auto_model consensus reads.fastq"  command? check if that exist and if it dose
// ============================================================
process MEDAKA_POLISH {
    label 'polish_env'
    tag   "medaka_polish"

    conda 'bioconda::medaka'
    container 'quay.io/biocontainers/medaka:2.2.2--py312h3050eb1_0'
    publishDir "${params.outdir}/04_polishing/medaka", mode: params.publish_mode

    input:
    path assembly
    path reads

    output:
    path "medaka_polished.fasta", emit: medaka_assembly
    path "versions.yml", emit: versions

    script:
    def explicit_model = params.medaka_model ?: ''
    """
    if [ -n "${explicit_model}" ]; then
        MODEL="${explicit_model}"
        echo "Using explicit params.medaka_model: \$MODEL"
    else
        # The basecaller model is in the read headers in one of two forms:
        #   - basecall_model_version_id=<model>        (MinKNOW / Guppy fastq)
        #   - RG:Z:<run-uuid>_<model>                  (Dorado SAM-style tags,
        #     e.g. RG:Z:741f927f-..._dna_r10.4.1_e8.2_400bps_sup@v5.2.0)
        # Scan the first 20000 reads (not just one) so a mix of models is caught.
        BCS=\$(gunzip -c ${reads} 2>/dev/null | awk 'NR%4==1' | head -n 20000 \\
              | tr ' \\t' '\\n\\n' \\
              | { grep -E '^(basecall_model_version_id=|RG:Z:)' || true; } \\
              | sed -E 's/^basecall_model_version_id=//; s/^RG:Z:[^_]*_//' \\
              | sort -u)
        NMOD=\$(printf '%s\\n' "\$BCS" | grep -c . || true)

        if [ "\$NMOD" -eq 0 ]; then
            echo "ERROR: no basecaller model (basecall_model_version_id= or RG:Z:) in the headers of ${reads}." >&2
            echo "       The medaka model materially affects consensus accuracy, so this" >&2
            echo "       pipeline will not pick one for you. Either supply reads with intact" >&2
            echo "       header metadata, or set --medaka_model explicitly." >&2
            echo "       Available names: medaka tools list_models" >&2
            exit 1
        elif [ "\$NMOD" -gt 1 ]; then
            echo "ERROR: reads carry \$NMOD different basecaller models:" >&2
            echo "\$BCS" >&2
            echo "       medaka needs exactly one; set --medaka_model explicitly." >&2
            exit 1
        fi
        BC="\$BCS"

        MODEL="\${BC}:consensus"
        echo "Detected basecaller model '\$BC' -> medaka model spec '\$MODEL'"
    fi

    # Validate BEFORE the expensive step. resolve_model errors out on an
    # unrecognised spec rather than quietly falling back to some default.
    medaka tools resolve_model --model "\$MODEL"

    medaka_consensus \\
      -i ${reads} \\
      -d ${assembly} \\
      -o medaka_out \\
      -t ${task.cpus} \\
      -m "\$MODEL"

    cp medaka_out/consensus.fasta medaka_polished.fasta
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: "\$(medaka --version 2>&1 | sed 's/medaka //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch medaka_polished.fasta
    echo "MEDAKA_POLISH stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: FUNANNOTATE_CLEAN
// Remove short contigs and small duplicated contigs.
// Clean and sort are separate processes so that a sort failure (e.g. a bad
// --base name) does not throw away the multi-hour clean result: the clean task
// stays cached and only the sort re-runs on -resume.
// ============================================================
process FUNANNOTATE_CLEAN {
    label 'polish_env'
    tag   "funannotate_clean"

    conda 'bioconda::funannotate'
    container 'quay.io/biocontainers/funannotate:1.8.17--pyhdfd78af_5'

    input:
    path assembly

    output:
    path "clean.fasta",  emit: clean_assembly
    path "versions.yml", emit: versions

    script:
    """
    # Clean: remove contigs shorter than 200 bp, and contigs shorter than N50
    # that are duplicated elsewhere (one minimap2 run per contig).
    # --cpus is the number of parallel minimap2 workers; funannotate's default is
    # 2, which is ~1 day for this assembly's ~16.6k small contigs. Each contig is
    # tested against a fixed reference, so the result does not depend on --cpus.
    funannotate clean -i ${assembly} -m 200 --cpus ${task.cpus} -o clean.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        funannotate: "\$(funannotate --version 2>&1 | sed 's/funannotate v//' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch clean.fasta
    echo "FUNANNOTATE_CLEAN stub completed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        funannotate: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: FUNANNOTATE_SORT
// Sort scaffolds by size (largest first) and rename them <base>_1..N
// ============================================================
process FUNANNOTATE_SORT {
    label 'polish_env'
    label 'res_light'
    tag   "funannotate_sort_${stage}"

    conda 'bioconda::funannotate'
    container 'quay.io/biocontainers/funannotate:1.8.17--pyhdfd78af_5'
    // Per-stage subdirectory: the polished assembly and (if the polish QC gate
    // falls back) the unpolished one are both published without clobbering.
    // (closure: the path refers to the `stage` input, so it must be evaluated per task)
    publishDir path: { "${params.outdir}/04_polishing/funannotate/${stage}" }, mode: params.publish_mode

    input:
    path assembly
    val  scaffold_basename
    val  stage      // 'polished' | 'unpolished'

    output:
    path "clean_sorted.fasta", emit: final_assembly
    path "versions.yml", emit: versions

    script:
    // funannotate sort refuses names over 16 characters (NCBI/GenBank limit),
    // i.e. <base>_<N> with N up to the scaffold count.
    // --minlen 0 (= no length filter, the documented default) must be explicit:
    // in 1.8.17 the argparse default is None and `None > 0` raises TypeError.
    """
    funannotate sort -i ${assembly} -o clean_sorted.fasta -b ${scaffold_basename} --minlen 0

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        funannotate: "\$(funannotate --version 2>&1 | sed 's/funannotate v//' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch clean_sorted.fasta
    echo "FUNANNOTATE_SORT stub completed"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        funannotate: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: BUSCO_FINAL
// Final BUSCO on polished assembly
// ============================================================
process BUSCO_FINAL {
    label 'polish_env'
    tag   "busco_${stage}_${lineage}"

    conda 'bioconda::busco'
    container 'quay.io/biocontainers/busco:6.1.0--pyhdfd78af_1'
    publishDir "${params.outdir}/04_polishing/qc/busco", mode: params.publish_mode

    input:
    path assembly
    val  lineage
    path busco_db_dir
    val  stage      // 'prepolish' | 'postpolish' -- keeps the two invocations' outputs distinct

    output:
    path "${stage}.${lineage}.full_table.tsv",     emit: busco_table
    path "${stage}.${lineage}.short_summary.txt",  emit: busco_summary
    path "${stage}.${lineage}.short_summary.json", emit: busco_summary_json
    path "versions.yml", emit: versions

    script:
    // BUSCO writes full_table.tsv under <outdir>/run_<lineage>/, and the short
    // summary as <outdir>/short_summary.specific.<lineage>.<outdir>.{txt,json}.
    // The previous version declared <outdir>/full_table.tsv and
    // <outdir>/short_summary.txt, neither of which BUSCO ever creates, so the
    // process always failed on missing outputs. Copy from the real paths, the
    // same way BUSCO_ASSEMBLY in modules/02_multi_assembler.nf does.
    //
    // --offline is required whenever an offline lineage dir is supplied: without
    // it BUSCO still contacts busco-data.ezlab.org to version-check and errors
    // out on a node with no egress, and a silent upstream dataset update would
    // make scores non-comparable across runs.
    //
    // _JAVA_OPTIONS: BUSCO runs bbtools stats.sh, whose default JVM heap is too
    // small for single-line (unwrapped) chromosome-scale FASTA records and dies
    // with OutOfMemoryError. Same fix as BUSCO_ASSEMBLY in 02_multi_assembler.nf.
    def busco_dir = params.busco_offline_dir ? "--offline --download_path ${busco_db_dir}" : ""
    def busco_mem = "${task.memory.toGiga()}g"
    """
    export _JAVA_OPTIONS="-Xmx${busco_mem}"
    busco -i ${assembly} -m geno -l ${lineage} -c ${task.cpus} \\
      -o busco_out ${busco_dir}

    cp busco_out/run_${lineage}/full_table.tsv                    ${stage}.${lineage}.full_table.tsv
    cp busco_out/short_summary.specific.${lineage}.busco_out.txt  ${stage}.${lineage}.short_summary.txt
    cp busco_out/short_summary.specific.${lineage}.busco_out.json ${stage}.${lineage}.short_summary.json
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: "\$(busco --version 2>&1 | sed 's/BUSCO //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    printf 'C:95.0%%[S:93.0%%,D:2.0%%],F:2.0%%,M:3.0%%,n:1000\\n' > ${stage}.${lineage}.short_summary.txt
    touch ${stage}.${lineage}.full_table.tsv ${stage}.${lineage}.short_summary.json
    echo "BUSCO_FINAL stub completed (${stage})"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: MERQURY_QV
// Reference-free consensus quality (QV) from the ONT read k-mer database.
//
// WHY: nothing in this pipeline previously checked that polishing made the
// assembly BETTER. BUSCO ran once, after the whole chain, so a polishing step
// that degraded the consensus was indistinguishable from one that helped.
// Merqury's QV compares assembly k-mers against the read k-mer spectrum, so it
// measures base-level accuracy directly and needs no reference genome.
//
// It reuses the meryl database stage 01 already builds -- no extra counting.
//
// CAVEAT for an ONT-only dataset: the read k-mer spectrum carries the platform's
// own error profile, so absolute QV is pessimistic compared with a
// short-read-derived spectrum. The DELTA between pre- and post-polish is the
// meaningful signal here, which is exactly what POLISH_QC_GATE compares.
// ============================================================
process MERQURY_QV {
    label 'polish_env'
    tag   "merqury_${stage}"

    conda 'bioconda::merqury=1.4.1'
    container 'quay.io/biocontainers/merqury:1.4.1--hdfd78af_1'
    publishDir "${params.outdir}/04_polishing/qc/merqury", mode: params.publish_mode

    input:
    path assembly
    path meryl_db
    val  stage      // 'prepolish' | 'postpolish'

    output:
    tuple val(stage), path("${stage}.qv"), emit: qv
    path "${stage}.completeness.stats", emit: completeness, optional: true
    path "${stage}_merqury_files/*",    emit: extras,       optional: true
    path "versions.yml", emit: versions

    script:
    """
    # merqury.sh does `source \$MERQURY/util/util.sh`, and MERQURY is only set by
    # conda activation, which does not happen under apptainer/docker exec. Unset,
    # it resolves to /util/util.sh and every helper (QV stats, plots) fails.
    export MERQURY=/usr/local/share/merqury

    # merqury.sh writes its outputs into the CWD prefixed with the 3rd argument.
    merqury.sh ${meryl_db} ${assembly} ${stage} || true

    # The .qv file is the only output the gate depends on, so it is the one thing
    # that must exist. Plot generation needs R and is allowed to fail without
    # taking the QV measurement down with it -- hence the `|| true` above.
    if [ ! -s ${stage}.qv ]; then
        echo "ERROR: merqury produced no QV for stage '${stage}'." >&2
        echo "       Check the merqury.sh output above (a missing /util/util.sh means MERQURY" >&2
        echo "       is unset), and that ${meryl_db} is a meryl >=1.4 database -- a 1.2-era" >&2
        echo "       database is not readable by meryl 1.4 and fails silently here." >&2
        exit 1
    fi

    mkdir -p ${stage}_merqury_files
    for f in ${stage}.*; do
        case "\$f" in
            ${stage}.qv|${stage}.completeness.stats) ;;
            *) mv "\$f" ${stage}_merqury_files/ 2>/dev/null || true ;;
        esac
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: "\$(meryl --version 2>&1 | sed 's/meryl //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    printf '${stage}\\t1000\\t500000\\t40.0000\\t1.0e-04\\n' > ${stage}.qv
    printf '${stage}\\tall\\t450000\\t500000\\t90.0000\\n'   > ${stage}.completeness.stats
    mkdir -p ${stage}_merqury_files
    touch ${stage}_merqury_files/.placeholder
    echo "MERQURY_QV stub completed (${stage})"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: POLISH_QC_GATE
// Compare pre- and post-polish QV and BUSCO completeness and decide which
// assembly continues downstream.
//
// A polishing step that degrades the consensus must not quietly propagate into
// the annotation and the paper. If either metric regresses beyond its tolerance
// the gate selects the UNPOLISHED assembly (polish_decision.txt = "unpolished")
// and the run continues with it; the polished assembly and this report are still
// published. Tolerances are params so a tiny, expected wobble is not a regression.
// A metric that cannot be parsed still FAILS the task: the gate refuses to pick
// an assembly it could not actually evaluate.
// ============================================================
process POLISH_QC_GATE {
    tag   "polish_qc_gate"

    conda 'conda-forge::python'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/04_polishing/qc", mode: params.publish_mode

    input:
    tuple val(pre_stage),  path(pre_qv)
    tuple val(post_stage), path(post_qv)
    path pre_busco_summary
    path post_busco_summary

    output:
    path "polish_qc_report.tsv", emit: report
    path "polish_decision.txt",  emit: decision   // "polished" | "unpolished"
    path "polish_qc_gate_mqc.tsv", emit: mqc      // MultiQC custom-content table
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import re, sys

    def read_qv(path):
        # merqury .qv columns: name, asm-only kmers, total kmers, QV, error rate
        for line in open(path):
            f = line.split()
            if len(f) >= 4:
                try:
                    return float(f[3])
                except ValueError:
                    continue
        return None

    def read_busco(path):
        m = re.search(r"C:([\\d.]+)%", open(path).read())
        return float(m.group(1)) if m else None

    pre_qv_v   = read_qv("${pre_qv}")
    post_qv_v  = read_qv("${post_qv}")
    pre_busco  = read_busco("${pre_busco_summary}")
    post_busco = read_busco("${post_busco_summary}")

    qv_tol    = float("${params.polish_qv_tolerance}")
    busco_tol = float("${params.polish_busco_tolerance}")

    rows, unparsed, regressions = [], [], []

    def check(label, pre, post, tol, unit):
        if pre is None or post is None:
            rows.append((label, "NA" if pre is None else pre,
                                "NA" if post is None else post, "", "UNPARSED"))
            unparsed.append(
                f"{label}: could not parse pre={pre} post={post} -- refusing to "
                f"pass a gate it cannot actually evaluate")
            return
        delta = post - pre
        ok = delta >= -tol
        rows.append((label, f"{pre:.4f}", f"{post:.4f}", f"{delta:+.4f}",
                     "PASS" if ok else "FAIL"))
        if not ok:
            regressions.append(
                f"{label} REGRESSED by {-delta:.4f} {unit} "
                f"({pre:.4f} -> {post:.4f}), beyond the {tol} {unit} tolerance")

    check("merqury_qv",         pre_qv_v,  post_qv_v,  qv_tol,    "QV")
    check("busco_completeness", pre_busco, post_busco, busco_tol, "pp")

    with open("polish_qc_report.tsv", "w") as out:
        out.write("metric\\tpre_polish\\tpost_polish\\tdelta\\tverdict\\n")
        for r in rows:
            out.write("\\t".join(str(x) for x in r) + "\\n")

    for r in rows:
        print("\\t".join(str(x) for x in r))

    if unparsed:
        sys.stderr.write("\\nERROR: polish QC gate cannot evaluate:\\n")
        for f in unparsed:
            sys.stderr.write("  - " + f + "\\n")
        sys.exit(1)

    if regressions:
        decision = "unpolished"
        sys.stderr.write("\\nWARNING: polishing degraded the assembly:\\n")
        for f in regressions:
            sys.stderr.write("  - " + f + "\\n")
        sys.stderr.write(
            "\\nContinuing with the UNPOLISHED assembly (it gets the same clean/sort).\\n"
            "The polished assembly is still published. To keep it instead, raise\\n"
            "--polish_qv_tolerance / --polish_busco_tolerance deliberately.\\n")
        print("\\npolish QC gate: REGRESSION -> using unpolished assembly")
    else:
        decision = "polished"
        print("\\npolish QC gate: PASS -> using polished assembly")

    with open("polish_decision.txt", "w") as d:
        d.write(decision + "\\n")

    # MultiQC custom content: the gate table, with the decision in the section
    # description so the report says plainly which assembly is the final product.
    final_busco = "postpolish" if decision == "polished" else "final_unpolished"
    desc = (f"Final assembly: {decision.upper()}. "
            f"Its BUSCO is the '{final_busco}' sample in the BUSCO section; "
            "'prepolish' and 'postpolish' are the gate's before/after measurements. "
            f"Tolerances: QV {qv_tol}, BUSCO {busco_tol} pp.")
    with open("polish_qc_gate_mqc.tsv", "w") as m:
        m.write("# id: 'polish_qc_gate'\\n")
        m.write("# section_name: '04 Polishing QC gate'\\n")
        m.write(f'# description: "{desc}"\\n')
        m.write("# plot_type: 'table'\\n")
        m.write("metric\\tpre_polish\\tpost_polish\\tdelta\\tverdict\\n")
        for r in rows:
            m.write("\\t".join(str(x) for x in r) + "\\n")

    # This script body is Python (see the shebang), so versions.yml is written
    # here -- a bash heredoc at this point is a Python SyntaxError.
    import platform
    with open("versions.yml", "w") as v:
        v.write('"${task.process}":\\n')
        v.write(f'    python: "{platform.python_version()}"\\n')
        v.write('    container: "${task.container ?: 'conda'}"\\n')
    """

    stub:
    """
    printf 'metric\\tpre_polish\\tpost_polish\\tdelta\\tverdict\\n'       > polish_qc_report.tsv
    printf 'merqury_qv\\t40.0000\\t42.0000\\t+2.0000\\tPASS\\n'         >> polish_qc_report.tsv
    printf 'busco_completeness\\t95.0000\\t95.4000\\t+0.4000\\tPASS\\n' >> polish_qc_report.tsv
    echo polished > polish_decision.txt
    cp polish_qc_report.tsv polish_qc_gate_mqc.tsv
    echo "POLISH_QC_GATE stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}
