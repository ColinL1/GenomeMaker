// ============================================================
// 09b_software_versions.nf
// Collate the per-process versions.yml files into the report's software table.
//
// Versions come from each tool's own --version at runtime, not container tags
// (tags can disagree with the tool they ship, and the conda profile has none).
// ============================================================

nextflow.enable.dsl = 2

process COLLATE_VERSIONS {
    label 'summary_env'
    tag   "collate_versions"

    // collectFile writes to a new path every run; 'deep' caches on content
    cache 'deep'

    conda 'conda-forge::python'
    container 'quay.io/biocontainers/python:3.12'
    publishDir "${params.outdir}/summary", mode: params.publish_mode

    input:
    path versions_file

    output:
    path "software_versions_mqc.yaml", emit: mqc
    path "software_versions.tsv",      emit: table
    path "versions.yml",               emit: versions

    script:
    """
    #!/usr/bin/env python3
    import glob, html, os, re

    # files written by this pipeline in a fixed two-level shape: no PyYAML needed
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
