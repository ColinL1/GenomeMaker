// ============================================================
// 09b_software_versions.nf
// Collate the per-process versions.yml files into the report's software table.
//
// WHY PER-PROCESS CAPTURE, not parsing container tags out of the module files:
//
//   1. A container TAG can disagree with the tool it ships. Verified on this
//      pipeline: `teambraker/braker3:v3.1.1` runs braker.pl **3.0.8**, and
//      `nanostat:0.1.5` runs NanoStat **0.0.1**. Parsing tags would have put the
//      wrong versions in the methods section.
//   2. Under `-profile conda` there are no containers at all, so a tag-based
//      approach reports nothing.
//
// Each process therefore runs its tool's own --version at runtime. Each command is
// wrapped in `|| echo unknown` so a bad command degrades to "unknown" rather than
// failing the task.
//
// The container is recorded ALONGSIDE the version, so the report carries both the
// tool's self-reported version and the exact image it came from.
// ============================================================

nextflow.enable.dsl = 2

process COLLATE_VERSIONS {
    label 'summary_env'
    tag   "collate_versions"

    // collectFile writes its output to a new work/tmp/ directory on every run, so
    // the default cache key (path + size + mtime) never matched and this task --
    // and MULTIQC_FINAL after it -- re-ran on every -resume. 'deep' keys on the
    // file's content instead.
    cache 'deep'

    conda 'conda-forge::python'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/summary", mode: params.publish_mode

    input:
    path versions_file      // one concatenated file; see subworkflows/09_summary.nf

    output:
    path "software_versions_mqc.yaml", emit: mqc
    path "software_versions.tsv",      emit: table
    path "versions.yml",               emit: versions

    script:
    """
    #!/usr/bin/env python3
    import glob, html, os, re

    # Minimal parser rather than a PyYAML dependency: the files this reads are
    # written by this pipeline in a fixed two-level shape.
    entries = {}
    proc, cur = None, {}
    for line in open("${versions_file}", errors="replace"):
        line = line.rstrip("\\n")
        if not line.strip():
            continue
        if not line.startswith(" ") and line.strip().endswith(":"):
            if proc and cur:
                entries.setdefault(proc, {}).update(cur)
            proc = line.strip().rstrip(":").strip('"')
            cur = {}
        elif ":" in line:
            k, v = line.split(":", 1)
            cur[k.strip()] = v.strip().strip('"')
    if proc and cur:
        entries.setdefault(proc, {}).update(cur)

    with open("software_versions.tsv", "w") as out:
        out.write("process\\ttool\\tversion\\tcontainer\\n")
        for proc in sorted(entries):
            e = entries[proc]
            cont = e.get("container", "")
            for tool, ver in sorted(e.items()):
                if tool == "container":
                    continue
                out.write(f"{proc}\\t{tool}\\t{ver}\\t{cont}\\n")

    # One row per distinct tool/version, with the images it came from.
    tools = {}
    for proc, e in entries.items():
        cont = e.get("container", "")
        for tool, ver in e.items():
            if tool == "container":
                continue
            tools.setdefault((tool, ver), set()).add(cont)

    rows = []
    for (tool, ver), conts in sorted(tools.items()):
        imgs = ", ".join(sorted(c for c in conts if c))
        rows.append(f"<tr><td><samp>{html.escape(tool)}</samp></td>"
                    f"<td>{html.escape(ver)}</td>"
                    f"<td><samp>{html.escape(imgs)}</samp></td></tr>")

    with open("software_versions_mqc.yaml", "w") as out:
        out.write("id: 'software_versions'\\n")
        out.write("section_name: 'Software versions'\\n")
        out.write("description: 'Each tool&apos;s OWN reported version, captured at runtime, "
                  "alongside the image it ran in. Container tags are not always accurate: in "
                  "this pipeline braker3:v3.1.1 reports braker.pl 3.0.8 and nanostat:0.1.5 "
                  "reports NanoStat 0.0.1.'\\n")
        out.write("plot_type: 'html'\\n")
        out.write("data: |\\n")
        out.write("    <table class='table'><thead><tr>"
                  "<th>Tool</th><th>Version</th><th>Container</th></tr></thead><tbody>\\n")
        for r in rows:
            out.write("    " + r + "\\n")
        out.write("    </tbody></table>\\n")

    # Self-report, so this process appears in its own table.
    with open("versions.yml", "w") as out:
        out.write('"COLLATE_VERSIONS":\\n')
        out.write(f'    python: "{".".join(str(x) for x in __import__("sys").version_info[:3])}"\\n')

    print(f"collated {len(entries)} processes, {len(tools)} distinct tool/version pairs")
    """

    stub:
    """
    printf 'process\\ttool\\tversion\\tcontainer\\n' > software_versions.tsv
    printf "id: 'software_versions'\\nsection_name: 'Software versions'\\nplot_type: 'html'\\ndata: |\\n    <p>stub</p>\\n" > software_versions_mqc.yaml
    printf '"COLLATE_VERSIONS":\\n    python: "stub"\\n' > versions.yml
    echo "COLLATE_VERSIONS stub completed"
    """
}
