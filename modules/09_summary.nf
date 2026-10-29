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

    # ---- polishing gate and the final assembly ----
    lines.append("")
    lines.append("POLISHING / FINAL ASSEMBLY")
    lines.append("-" * 68)
    if os.path.exists("polish_decision.txt"):
        decision = open("polish_decision.txt").read().strip()
        lines.append(f"  Polish QC gate: {decision} contigs carried forward")
    else:
        lines.append("  polish QC gate did not run")
    for stage, label in (("prepolish", "before polishing (gate input)"),
                         ("postpolish", "after polishing (gate input)"),
                         ("final", "final assembly")):
        for f in sorted(glob.glob(f"{stage}.*short_summary*.txt")):
            m = re.search(r"C:([\\d.]+)%\\[S:([\\d.]+)%,D:([\\d.]+)%\\],F:([\\d.]+)%,M:([\\d.]+)%",
                          open(f, errors="replace").read())
            if m:
                lines.append(f"  BUSCO {label:<30} C={m.group(1)}% D={m.group(3)}% M={m.group(5)}%")

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
    def config_arg = (multiqc_config && multiqc_config.baseName != 'NO_FILE') ? "-c ${multiqc_config}" : ""
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

// Optional file input: '--flag path', or nothing for a NO_ placeholder
def opt_flag(flag, f) {
    return (f && !f.baseName.startsWith('NO_')) ? "${flag} ${f}" : ''
}

// ============================================================
// Process: GENOME_SUMMARY
// summary.tsv (bin/genome_summary.py): assembly, completeness, annotation, repeats and
// mitogenome numbers, computed from the delivered files. Run on final/ (strategy 1) and, with
// the release curation, on release_curation/ (that one is the MultiQC summary table then).
// ============================================================
process GENOME_SUMMARY {
    label 'res_light'
    tag   "genome_summary_${label}"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir path: { "${params.outdir}/summary${params.release_folder_suffix}/${label}" }, mode: params.publish_mode

    input:
    val  label                                   // 'final' | 'release_curation'
    path fasta,          stageAs: 'genome/*'
    path anno,           stageAs: 'anno/*'       // funannotate output dir, or NO_ placeholder
    path gff,            stageAs: 'gff/*'        // annotation GFF3 to count from (or NO_: anno's GFF3)
    path compleasm,      stageAs: 'compleasm/*'
    path busco,          stageAs: 'busco/*'
    path qv,             stageAs: 'qv/*'
    path repeat_tbl,     stageAs: 'repeats/*'
    path hc_report,      stageAs: 'hc/*'
    path mito_fasta,     stageAs: 'mito/*'
    path mito_check,     stageAs: 'mitocheck/*'
    path symbiont_busco, stageAs: 'symbiont/*'

    output:
    path "summary.tsv", emit: summary

    script:
    def has_anno = anno && !anno.baseName.startsWith('NO_')
    def gff_arg  = (gff && !gff.baseName.startsWith('NO_')) ? "--gff ${gff}" : (has_anno ? '--gff "$(ls anno/*/annotate_results/*.gff3 | grep -v \'\\.hc\\.\' | head -1)"' : '')
    def anno_args = has_anno
        ? '--stats "$(ls anno/*/annotate_results/*.stats.json | head -1)" --annotations "$(ls anno/*/annotate_results/*.annotations.txt | head -1)"'
        : ''
    """
    genome_summary.py --fasta ${fasta} ${gff_arg} ${anno_args} \\
      ${opt_flag('--compleasm', compleasm)} ${opt_flag('--busco', busco)} ${opt_flag('--qv', qv)} \\
      ${opt_flag('--repeat-tbl', repeat_tbl)} ${opt_flag('--hc-report', hc_report)} \\
      ${opt_flag('--mito-fasta', mito_fasta)} ${opt_flag('--mito-check', mito_check)} \\
      ${opt_flag('--symbiont-busco', symbiont_busco)} --funannotate-busco-db '${params.funannotate_busco_db}' \\
      --n-chrom ${params.n_chromosome_scaffolds ?: 0} --out summary.tsv
    """

    stub:
    """
    printf 'section\\tmetric\\tvalue\\nAssembly\\tTotal length (bp)\\t1\\nAnnotation\\tProtein-coding genes\\t1\\n' > summary.tsv
    """
}

// ============================================================
// Process: METHODS
// The methods text of the run (bin/methods.py): the template main.nf builds from the steps
// that ran, filled with the recorded tool versions and the automatic decisions read from
// the decision files. Shown at the bottom of the MultiQC report and written to qc/methods.md.
// ============================================================
process METHODS {
    label 'res_light'
    tag   "methods"

    conda 'conda-forge::python=3.12'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/summary${params.release_folder_suffix}", mode: params.publish_mode

    input:
    path template
    path versions
    path ranking,         stageAs: 'd/ranking/*'
    path polish_report,   stageAs: 'd/polish_report/*'
    path polish_decision, stageAs: 'd/polish_decision/*'
    path purge_report,    stageAs: 'd/purge/*'
    path mito_removed,    stageAs: 'd/mito/*'
    path length_removed,  stageAs: 'd/length/*'
    path screen_actions,  stageAs: 'd/screen/*'

    output:
    path "methods.md", emit: methods

    script:
    """
    methods.py ${template} ${versions} methods.md --ranking ${ranking} --polish-report ${polish_report} \\
      --polish-decision ${polish_decision} --purge-report ${purge_report} --mito-removed ${mito_removed} \\
      --length-removed ${length_removed} --screen-actions ${screen_actions}
    """

    stub:
    """
    cp ${template} methods.md
    """
}

// ============================================================
// Process: FINAL_ASSEMBLE
// <outdir>/final/: the automatic product, written in both strategies:
//   genome/      final assembly (soft-masked when repeats ran), scaffold_name_map.tsv (with RagTag)
//   annotation/  funannotate GFF3, proteins, CDS and mRNA transcripts, annotations.txt (stage 08)
//   repeats/     repeats.gff3, RepeatMasker table, classified library
//   mitogenome/  sequence and annotation (stage 05)
//   symbiont/    the symbiont bin of the BlobToolKit split
//   qc/          summary.tsv, methods.md and the decision files of the automatic gates
//   multiqc_report.html
// ============================================================
process FINAL_ASSEMBLE {
    label 'res_light'
    tag   "final_assemble"

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path fasta,         stageAs: 'in/genome/assembly.fa'
    path name_map,      stageAs: 'in/genome/map/*'
    path anno,          stageAs: 'in/anno/*'
    path repeat_files,  stageAs: 'in/repeats/*'
    path mito_files,    stageAs: 'in/mito/*'
    path symbiont,      stageAs: 'in/symbiont/*'
    path qc_files,      stageAs: 'in/qc/*'
    path report,        stageAs: 'in/report/*'

    output:
    path "final", emit: dir

    script:
    """
    set -euo pipefail
    mkdir -p final/{genome,qc}
    cp -L in/genome/assembly.fa final/genome/${params.species_name}.fa
    for f in in/genome/map/* in/qc/* in/repeats/* in/mito/* in/symbiont/*; do
        [ -e "\$f" ] || continue
        case \$(basename \$f) in NO_*) continue ;; esac
        case \$f in
            in/genome/map/*) cp -L \$f final/genome/ ;;
            in/qc/*)         cp -L \$f final/qc/ ;;
            in/repeats/*)    mkdir -p final/repeats    && cp -L \$f final/repeats/ ;;
            in/mito/*)       mkdir -p final/mitogenome && cp -L \$f final/mitogenome/ ;;
            in/symbiont/*)   mkdir -p final/symbiont   && cp -L \$f final/symbiont/symbiont_bin.fa ;;
        esac
    done
    for d in in/anno/*; do
        [ -d "\$d/annotate_results" ] || continue
        mkdir -p final/annotation
        for f in \$d/annotate_results/*.{gff3,proteins.fa,cds-transcripts.fa,mrna-transcripts.fa,annotations.txt}; do
            if [ -e "\$f" ]; then cp -L \$f final/annotation/; fi
        done
    done
    cp -L in/report/multiqc_report.html final/
    """
    // no stub: the script only copies files, so -stub runs check the folder layout
}

// ============================================================
// Process: PUBLISH_MULTIQC
// Copies the MultiQC report into the release-curation folder (<outdir>/release_curation/<name>/),
// which is built before MultiQC runs.
// ============================================================
process PUBLISH_MULTIQC {
    label 'res_light'
    tag   "publish_multiqc"

    publishDir "${params.outdir}/release_curation${params.release_folder_suffix}/${params.release_name}", mode: 'copy'

    input:
    path report, stageAs: 'in/*'
    path release_dir         // the release-curation folder: copied only after it was built

    output:
    path "multiqc_report.html"

    script:
    """
    cp -L in/multiqc_report.html multiqc_report.html
    """
}
