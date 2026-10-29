// ============================================================
// 09_summary.nf
// Aggregate all QC results into a final summary report
// ============================================================

nextflow.enable.dsl = 2

// ============================================================
// Process: GENERATE_SUMMARY
// Collect gfastats, BUSCO, QUAST, annotation stats into reports
// ============================================================
process GENERATE_SUMMARY {
    label 'res_small_2h'
    tag   "generate_summary"

    // Plain Python: the script uses only the standard library. The previous
    // mulled pandas/jinja2 image is an R environment with no python3 at all
    // ("env: can't execute 'python3'").
    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'

    publishDir "${params.outdir}/summary${params.release_folder_suffix}", mode: params.publish_mode

    input:
    path gfastats_files
    path busco_files
    path polish_files       // polishing-stage BUSCO short summaries + polish_decision.txt
    path comparison_tsv
    path comparison_json
    path annotation_stats

    output:
    path "summary/assembly_comparison.tsv", emit: summary_tsv
    path "summary/assembly_comparison.json", emit: summary_json
    path "summary/pipeline_report.txt",      emit: pipeline_report
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    # Produces a data-bearing plain-text summary.
    #
    # The previous version wrote a hardcoded heredoc whose every section said
    # "(See assembly_comparison.tsv for details)" -- it contained no data at all,
    # and was indented, so the output was indented too. It looked like a report
    # without being one.
    #
    # MultiQC is the full interactive report; this is the quick terminal/methods
    # read, generated from the same JSON MultiQC consumes.
    import glob, json, os, re

    os.makedirs("summary", exist_ok=True)

    def first(pat):
        hits = sorted(glob.glob(pat))
        return hits[0] if hits else None

    comp = first("*comparison*.json")
    stats = first("*annotation_stats*.json") or first("*stats*.json")

    lines = []
    lines.append("=" * 68)
    lines.append("Coral Genome Assembly & Annotation -- Summary")
    lines.append("=" * 68)

    # ---- per-assembler comparison ----
    lines.append("")
    lines.append("ASSEMBLY COMPARISON")
    lines.append("-" * 68)
    if comp:
        data = json.load(open(comp))
        lines.append(f"{'assembler':<16}{'BUSCO C%':>9}{'dup%':>7}{'N50':>16}{'total bp':>16}")
        for name in sorted(data):
            b = (data[name] or {}).get("busco", {}) or {}
            g = (data[name] or {}).get("gfastats", {}) or {}
            def gf(*keys):
                for k in keys:
                    for gk, gv in g.items():
                        if gk.lower() == k.lower():
                            return gv
                return "-"
            lines.append(f"{name:<16}{b.get('complete_pct','-'):>9}{b.get('duplicated_pct','-'):>7}"
                         f"{gf('Scaffold N50','Contig N50'):>16}{gf('Total scaffold length','Total contig length'):>16}")
    else:
        lines.append("  no comparison JSON was staged")

    # ---- polishing: which assembly is the final product ----
    lines.append("")
    lines.append("POLISHING / FINAL ASSEMBLY")
    lines.append("-" * 68)
    decision = None
    if os.path.exists("polish_decision.txt"):
        decision = open("polish_decision.txt").read().strip()
        lines.append(f"  Final assembly (polish QC gate): {decision.upper()}")
    else:
        lines.append("  polish QC gate did not run: the polished assembly is the final product")
    final_stage = "final_unpolished" if decision == "unpolished" else "postpolish"
    for stage, label in (("prepolish", "before polishing (gate input)"),
                         ("postpolish", "after polishing (gate input)"),
                         ("final_unpolished", "unpolished, cleaned+sorted")):
        for f in sorted(glob.glob(f"{stage}.*short_summary*.txt")):
            m = re.search(r"C:([\\d.]+)%\\[S:([\\d.]+)%,D:([\\d.]+)%\\],F:([\\d.]+)%,M:([\\d.]+)%",
                          open(f, errors="replace").read())
            if m:
                tag = "   <-- FINAL PRODUCT" if stage == final_stage else ""
                lines.append(f"  BUSCO {label:<30} C={m.group(1)}% D={m.group(3)}% M={m.group(5)}%{tag}")

    # ---- BUSCO summaries actually present ----
    lines.append("")
    lines.append("BUSCO SUMMARIES FOUND")
    lines.append("-" * 68)
    found = False
    for f in sorted(glob.glob("*short_summary*.txt")) + sorted(glob.glob("*full_table*.tsv")):
        m = re.search(r"C:([\\d.]+)%\\[S:([\\d.]+)%,D:([\\d.]+)%\\],F:([\\d.]+)%,M:([\\d.]+)%",
                      open(f, errors="replace").read())
        if m:
            found = True
            lines.append(f"  {f:<46} C={m.group(1)}% D={m.group(3)}% M={m.group(5)}%")
    if not found:
        lines.append("  none parseable")

    # ---- annotation ----
    lines.append("")
    lines.append("ANNOTATION")
    lines.append("-" * 68)
    if stats:
        try:
            d = json.load(open(stats))
            flat = d.get("annotation", d) if isinstance(d, dict) else {}
            for k in sorted(flat):
                v = flat[k]
                if not isinstance(v, (dict, list)):
                    lines.append(f"  {k:<40} {v}")
        except Exception as e:
            lines.append(f"  could not parse {stats}: {e}")
    else:
        lines.append("  functional annotation did not run, or produced no stats JSON")

    lines.append("")
    lines.append("Full interactive report: summary/multiqc/multiqc_report.html")
    lines.append("Software provenance:     summary/software_versions.tsv")
    lines.append("=" * 68)

    open("summary/pipeline_report.txt", "w").write("\\n".join(lines) + "\\n")

    # Carry the comparison tables through unchanged.
    import shutil
    for src, dst in (("${comparison_tsv}", "summary/assembly_comparison.tsv"),
                     ("${comparison_json}", "summary/assembly_comparison.json")):
        try:
            shutil.copyfile(src, dst)
        except Exception:
            pass

    print("\\n".join(lines))

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
    mkdir -p summary
    touch summary/assembly_comparison.tsv summary/assembly_comparison.json summary/pipeline_report.txt
    echo "GENERATE_SUMMARY stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: MULTIQC_CUSTOM_CONTENT
// Formats the outputs MultiQC has no module for (gfastats, Merqury, blobtools,
// RepeatMasker, barrnap, BRAKER4 QC, tRNAscan, Infernal, StringTie, compleasm,
// OMArk, GFF3 validation, funannotate, Phobius) as MultiQC custom content
// (*_mqc.tsv / *_mqc.png). Section ids start with their stage (stageNN_) so
// multiqc_config.yaml can order the report 00 -> 09. Script: bin/.
// ============================================================
process MULTIQC_CUSTOM_CONTENT {
    label 'res_small'
    tag   "multiqc_custom_content"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'

    input:
    path report_files, stageAs: 'inputs/*'

    output:
    path "mqc/*", emit: mqc
    path "versions.yml", emit: versions

    script:
    """
    # gene models are refined before OMArk/StringTie in this pipeline (stage 07c)
    MQC_REFINED_UPSTREAM=1 multiqc_custom_content.py inputs mqc

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "\$(python3 --version 2>&1 | sed 's/Python //')"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p mqc
    printf '# id: stub\\nA\\tB\\n' > mqc/stub_mqc.tsv
    printf 'Metric\\tValue\\nProtein-coding genes\\t1\\n' > mqc/stage08_funannotate_mqc.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
    END_VERSIONS
    """
}

// ============================================================
// Process: MULTIQC
// Aggregate every QC-producing step of the pipeline into a
// single MultiQC report: read QC (pre-assembly), per-assembler
// gfastats/QUAST/BUSCO, the post-polish BUSCO re-run, and
// functional annotation stats (e.g. funannotate/InterProScan/
// eggNOG summary json, if MultiQC-recognisable).
//
// This process is intentionally permissive about its inputs:
// it just stages whatever files it is given into one directory
// and lets `multiqc` auto-detect what it recognises. Missing or
// empty optional channels (e.g. mitogenome/annotation not run)
// are fine — MultiQC will simply report on what's present.
// ============================================================
process MULTIQC_FINAL {
    label 'res_tiny'
    tag   "multiqc"

    conda 'bioconda::multiqc'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'

    publishDir "${params.outdir}/summary${params.release_folder_suffix}/multiqc", mode: params.publish_mode

    input:
    path qc_files, stageAs: 'qc_inputs/*'
    path multiqc_config

    output:
    path "multiqc_report.html", emit: report
    // MultiQC names the data dir after --filename: multiqc_report.html -> multiqc_report_data
    path "multiqc_report_data", emit: data
    path "multiqc_report_data.zip", emit: data_zip, optional: true
    path "versions.yml", emit: versions

    script:
    def config_arg = (multiqc_config && multiqc_config.name != 'NO_FILE') ? "-c ${multiqc_config}" : ""
    """
    multiqc qc_inputs \\
        --force \\
        --filename multiqc_report.html \\
        --outdir . \\
        ${config_arg}

    # Keep a zipped copy of the parsed data alongside the html report
    if [ -d multiqc_report_data ]; then
        # zip is not in the MultiQC image; Python (which it has) writes the same archive.
        python3 -m zipfile -c multiqc_report_data.zip multiqc_report_data/
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: "\$(multiqc --version 2>&1 | sed 's/.*version //' || echo unknown)"
        container: "${task.container ?: 'conda'}"
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html
    mkdir -p multiqc_report_data
    echo stub > multiqc_report_data/stub.txt
    touch multiqc_report_data.zip
    echo "MULTIQC stub completed"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: "stub"
    END_VERSIONS
    """
}
