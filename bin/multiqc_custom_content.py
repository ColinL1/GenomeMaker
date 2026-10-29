#!/usr/bin/env python3
"""Convert pipeline outputs that MultiQC has no module for into MultiQC
custom-content files (*_mqc.tsv tables, *_mqc.png images).

Usage: multiqc_custom_content.py <input_dir> <output_dir>

Every section id starts with its pipeline stage (stage02_..., stage07_...) so
multiqc_config.yaml can order the report by stage. A missing input simply means
that section is not written (e.g. a stage that was disabled). Standard library only.
"""
import base64
import glob
import gzip
import json
import os
import re
import shutil
import sys
from collections import Counter, defaultdict

IN, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)


def find(pattern):
    """All files under IN matching a glob (recursive, sorted)."""
    return sorted(glob.glob(os.path.join(IN, "**", pattern), recursive=True))


def table(sec_id, name, description, header, rows, fmt=None,
          parent_id=None, parent_name=None, parent_description=None):
    """Write a MultiQC custom-content table (first column = row name)."""
    if not rows:
        return
    with open(os.path.join(OUT, f"{sec_id}_mqc.tsv"), "w") as fh:
        fh.write(f"# id: '{sec_id}'\n")
        fh.write(f"# section_name: '{name}'\n")
        fh.write(f'# description: "{description}"\n')
        if parent_id:
            fh.write(f"# parent_id: '{parent_id}'\n")
            fh.write(f"# parent_name: '{parent_name}'\n")
            fh.write(f'# parent_description: "{parent_description}"\n')
        fh.write("# plot_type: 'table'\n")
        fh.write("# pconfig:\n")
        fh.write(f"#     id: '{sec_id}_table'\n")
        fh.write("#     col1_header: '" + header[0] + "'\n")
        if fmt:
            fh.write(f"#     format: '{fmt}'\n")
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def image(sec_id, src):
    """Copy a PNG as <id>_mqc.png; its section name comes from multiqc_config."""
    if src and os.path.exists(src):
        shutil.copyfile(src, os.path.join(OUT, f"{sec_id}_mqc.png"))


def plot(sec_id, name, description, plot_type, pconfig, data, cats=None,
         parent_id=None, parent_name=None, parent_description=None):
    """Write a MultiQC custom-content JSON file (linegraph/bargraph/scatter/violin).

    `data` is the plot_type-specific MultiQC custom-content shape, e.g. for
    linegraph {series: {x: y}}, for bargraph {sample: {category: value}}.
    """
    if not data:
        return
    payload = {
        "id": sec_id,
        "section_name": name,
        "description": description,
        "plot_type": plot_type,
        "pconfig": {"id": f"{sec_id}_plot", **pconfig},
        "data": data,
    }
    if cats is not None:
        payload["cats"] = cats
    if parent_id:
        payload["parent_id"] = parent_id
        payload["parent_name"] = parent_name
        payload["parent_description"] = parent_description
    with open(os.path.join(OUT, f"{sec_id}_mqc.json"), "w") as fh:
        json.dump(payload, fh, indent=1)


def html(sec_id, name, description, body, parent_id=None, parent_name=None, parent_description=None):
    """Write a MultiQC custom-content section of raw HTML (e.g. hand-laid-out images)."""
    payload = {"id": sec_id, "section_name": name, "description": description,
               "plot_type": "html", "data": body}
    if parent_id:
        payload["parent_id"] = parent_id
        payload["parent_name"] = parent_name
        payload["parent_description"] = parent_description
    with open(os.path.join(OUT, f"{sec_id}_mqc.json"), "w") as fh:
        json.dump(payload, fh, indent=1)


def img_tag(src, min_width="480px", max_width="700px"):
    """A responsive flex-item <img> with the PNG inlined as base64 (same
    mechanism MultiQC uses for bare *_mqc.png files). Meant to sit in a
    `display:flex; flex-wrap:wrap` container: items this wide side by side
    when there's room, wrapping to one per row (full width, up to max_width)
    when the viewport is too narrow for both - no media query needed."""
    b64 = base64.b64encode(open(src, "rb").read()).decode("ascii")
    item_style = f"flex:1 1 {min_width}; min-width:{min_width}; max-width:{max_width};"
    return (f'<div style="{item_style}"><img src="data:image/png;base64,{b64}" '
            f'style="width:100%; height:auto; display:block;" /></div>')


def decimate(points, keep_dense=200, stride=20):
    """Thin a long, smooth (x, y) curve for plotting: keep every point in the
    fast-changing head, then every `stride`-th point, always keeping the last."""
    if len(points) <= keep_dense:
        return points
    out = points[:keep_dense] + points[keep_dense::stride]
    if out[-1] != points[-1]:
        out.append(points[-1])
    return out


def final_gff():
    """The released annotation: funannotate's annotate_results GFF3, else the refined GFF3."""
    fin = [f for f in find("*.gff3") if "annotate_results" in f and ".hc." not in f]
    return (fin or find("refined.gff3") or [None])[0]


def refined_gff():
    """The refined (pre-funannotate) GFF3; it keeps ncRNA_class, which funannotate drops."""
    return (find("refined.gff3") or [None])[0]


def protein_coding(a):
    """Protein-coding genes from funannotate's stats.json 'annotation' block."""
    g = a.get("genes")
    if g is None:
        return None
    pc = g - a.get("tRNA", 0) - a.get("rRNA", 0) - a.get("ncRNA", 0)
    if pc <= 0:
        raise SystemExit(f"protein-coding gene count {pc} <= 0 from stats.json")
    return pc


def stem(path, suffix):
    return os.path.basename(path)[: -len(suffix)] if path.endswith(suffix) else os.path.basename(path)


# ---- 02 assembly contiguity (gfastats, every assembly stage) --------------------
GFA_ORDER = ["input_assembly", "purged", "filtered", "prepolish", "postpolish"]
rows = []
for f in find("*.gfastats.txt"):
    kv = {}
    for line in open(f):
        if ":" in line:
            k, v = line.split(":", 1)
            kv[k.strip()] = v.strip()
    s = stem(f, ".gfastats.txt")
    rows.append((s, kv.get("# scaffolds", "NA"), kv.get("Total scaffold length", "NA"),
                 kv.get("Scaffold N50", "NA"), kv.get("# contigs", "NA"),
                 kv.get("Contig N50", "NA"), kv.get("# gaps in scaffolds", "NA"),
                 kv.get("GC content %", "NA")))
rows.sort(key=lambda r: GFA_ORDER.index(r[0]) if r[0] in GFA_ORDER else 99)
table("stage02_gfastats", "02 Assembly: contiguity (gfastats)",
      "Assembly statistics at every stage: input, after purge_dups, after the taxonomic filter, before and after polishing.",
      ["Stage", "Scaffolds", "Total length (bp)", "Scaffold N50", "Contigs", "Contig N50", "Gaps", "GC %"], rows)

# ---- 03 taxonomic filter (blobtools) -------------------------------------------
rows = []
for f in find("*_summary.json"):
    cat = stem(f, "_summary.json")
    if cat not in ("host", "symbiont"):
        continue
    try:
        tot = json.load(open(f))["summaryStats"]["hits"]["total"]
    except Exception:
        continue
    gc, cov = tot.get("gc", ["NA"]), tot.get("cov", ["NA"])
    rows.append((cat, tot.get("count", "NA"), tot.get("span", "NA"), tot.get("n50", "NA"),
                 round(gc[0] * 100, 2) if isinstance(gc[0], (int, float)) else "NA",
                 cov[0] if cov else "NA"))
PARENT_FILTERING = dict(
    parent_id="stage03_filtering", parent_name="03 Filtering",
    parent_description="Taxonomic filtering with BlobToolKit.")
table("stage03_blobtools_filter", "Taxonomic filter: category table",
      "Scaffolds kept per category by the blobtools filter (host = Cnidaria + no-hit; symbiont = Symbiodiniaceae).",
      ["Category", "Scaffolds", "Span (bp)", "N50", "GC %", "Mean coverage"], rows,
      **PARENT_FILTERING)

# blob + snail plots: still static PNGs (BlobToolKit's own renderer); per-scaffold
# GC/coverage isn't exported by this pipeline stage, so they can't be re-plotted
# natively. Side by side (each min 480px, so still readable) when there's room;
# flex-wrap drops the second one to its own full-width row otherwise.
imgs = []
for kind in ("blob.circle", "snail"):
    hit = find(f"*.{kind}.png")
    if hit:
        imgs.append(img_tag(hit[0]))
if imgs:
    body = ('<div style="display:flex; flex-wrap:wrap; gap:16px; justify-content:center;">'
            + "".join(imgs) + "</div>")
    html("stage03_blobtools_plots", "BlobToolKit blob + snail plots",
         "Blob plot (GC vs coverage, bubble size = scaffold span, coloured by taxonomy) and "
         "snail plot (contiguity, composition, BUSCO), host-filtered scaffolds.",
         body, **PARENT_FILTERING)

# Cumulative span: native linegraph from the scaffold id+length lists (host vs
# symbiont — the finer Cnidaria/no-hit split used by the original PNG is not
# recoverable without per-scaffold taxonomy, which isn't exported either).
def cumulative_series(path):
    lens = []
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        next(fh, None)  # header
        for line in fh:
            c = line.rstrip("\r\n").split("\t")
            if len(c) >= 3 and c[2].isdigit():
                lens.append(int(c[2]))
    lens.sort(reverse=True)
    points, cum = [], 0
    for i, l in enumerate(lens, 1):
        cum += l
        points.append((i, cum))
    return decimate(points)

# NOTE: pass [x, y] pairs, not a {x: y} dict. MultiQC's custom-content loader
# tries to coerce dict-string-keys to float for sorting, but (at least in
# MultiQC 1.33) the coerced result is assigned to the wrong variable and never
# reaches the renderer, so string keys get sorted lexicographically ("10000"
# before "2000") -> a zigzagging line instead of a smooth curve. Pairs avoid
# that path entirely since the numbers are never stringified to begin with.
cum_data = {}
for cat, fname in (("host", "host_scaffolds.tsv"), ("symbiont", "symbiont_scaffolds.tsv")):
    hit = find(fname)
    series = cumulative_series(hit[0]) if hit else []
    if series:
        cum_data[cat] = series
plot("stage03_blobtools_cumulative", "Cumulative span", "Cumulative scaffold span by rank "
     "(longest first), host vs symbiont.", "linegraph",
     {"xlab": "scaffold rank (longest first)", "ylab": "cumulative length (bp)", "x_decimals": 0},
     cum_data, **PARENT_FILTERING)

# ---- 03c contamination screens (FCS-adaptor, VecScreen, FCS-GX) ------------------
# screen_actions.tsv (apply_screen_actions.py) and screen_review.tsv (screen_review.py).
for f in find("screen_actions.tsv")[:1]:
    cnt = Counter()
    for line in open(f):
        x = line.rstrip("\n").split("\t")
        if x[0] in ("applied", "not_applied"):
            cnt[(x[0], x[2])] += 1
    rows = [(f"{a} {s}", n) for (s, a), n in sorted(cnt.items())]
    table("stage03c_screen_actions", "Contamination screens: actions",
          "Actions from FCS-GX, FCS-adaptor and the review, applied to the contigs before scaffolding "
          "(EXCLUDE drops a sequence, TRIM/FIX remove or mask a range).", ["Action", "Count"], rows)
for f in find("screen_review.tsv")[:1]:
    cnt = Counter()
    for line in open(f):
        x = line.rstrip("\n").split("\t")
        if x[0] != "seq" and len(x) > 6:
            cnt[(x[3], x[4], x[6])] += 1
    rows = [(f"{src} {cat}: {d}", n) for (src, cat, d), n in sorted(cnt.items())]
    table("stage03c_screen_review", "Contamination screens: reviewed hits",
          "FCS-GX REVIEW and VecScreen Strong/Moderate hits, with the decision taken for each.",
          ["Hit", "Count"], rows)

# ---- 04 polishing: Merqury ------------------------------------------------------
qv, comp = {}, {}
for f in find("*.qv"):
    parts = open(f).readline().split()
    if len(parts) >= 5:
        qv[stem(f, ".qv")] = (parts[3], parts[4])
for f in find("*.completeness.stats"):
    parts = open(f).readline().split()
    if len(parts) >= 5:
        comp[stem(f, ".completeness.stats")] = parts[4]
rows = [(s, qv[s][0], qv[s][1], comp.get(s, "NA")) for s in sorted(qv, key=lambda x: (x != "prepolish", x))]
PARENT_POLISHING = dict(
    parent_id="stage04_polishing", parent_name="04 Polishing",
    parent_description="Consensus polishing (medaka) and its QC gate.")
table("stage04_merqury", "Merqury QV and k-mer completeness",
      "Consensus quality (QV) and k-mer completeness of the assembly before and after polishing (read k-mers: Q15 ONT).",
      ["Stage", "QV", "Error rate", "k-mer completeness %"], rows, **PARENT_POLISHING)
# Note: the polish QC gate table (polish_qc_gate_mqc.tsv) is written by an earlier
# pipeline step, not this script, so its parent_id can't be set from here. To fully
# group "04 Polishing", add parent_id: 'stage04_polishing' (matching the id/name/
# description above) to that step's own mqc-table header.

# ---- 05 mitogenome: sequence source, protein check, liftover + MITOS2 annotation ----
PARENT_MITO = dict(
    parent_id="stage05_mitogenome", parent_name="05 Mitogenome",
    parent_description="Mitochondrial genome: sequence, annotation by liftover and by MITOS2. "
                       "Separate from the nuclear assembly.")


def kv_file(pattern):
    for f in find(pattern):
        return dict(l.rstrip("\n").split("\t", 1) for l in open(f) if "\t" in l)
    return {}


def gff_rows(path):
    """(name, type, start, end, strand) per gene-level feature of a GFF3."""
    feats = []
    for line in open(path):
        c = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(c) < 9:
            continue
        feats.append((c[2], int(c[3]), int(c[4]), c[6],
                      dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)))
    return feats


mito_summary = []   # reused by the pipeline summary
src = kv_file("mitogenome_source.tsv")
check_rows = []
for f in find("mito_gene_check.tsv"):
    check_rows = [l.rstrip("\n").split("\t") for l in open(f)][1:]
pol = kv_file("mito_polish_stats.tsv")
getorg_logs = find("get_org.log.txt")
if src or check_rows or pol or getorg_logs:
    # without an annotation reference there is no rotation step, so no source file
    source = src.get("source") or ("Illumina DNA, GetOrganelle" if getorg_logs
                                   else "ONT contig (MitoHiFi), Medaka-polished")
    status = re.search(r"Result status of animal_mt: (.+)", open(getorg_logs[0]).read()) if getorg_logs else None
    no_mito = bool(getorg_logs) and not (status and status.group(1).strip() == "circular genome")
    if no_mito:
        source += " -- NO complete circular mitogenome was produced; nothing was annotated"
    mito_summary = [("Sequence source", source),
                    ("Length (bp)", src.get("length") or pol.get("polished_length", "NA")),
                    ("Annotation reference (liftover, rotation origin)",
                     src.get("annotation_reference", "not run" if no_mito else "none (MITOS2 only)"))]
    for f in getorg_logs:
        txt = open(f).read()
        for label, pat in (("GetOrganelle result", r"Result status of animal_mt: (.+)"),
                           ("GetOrganelle read pairs used", r"Reads used = (\S+)"),
                           ("GetOrganelle mitogenome base coverage", r"Average animal_mt base-coverage = (\S+)"),
                           ("Illumina input files", r"concatenated in this order \(R1 then R2 lists\): (.+)")):
            m = re.search(pat, txt)
            if m:
                mito_summary.append((label, m.group(1).strip()))
    if pol:
        depth = int(pol.get("bases_used", 0)) / max(int(pol.get("polished_length", 1)), 1)
        mito_summary += [("Medaka polish: ONT mito reads used",
                          f'{pol.get("reads_used", "NA")} reads of mean quality >= Q{pol.get("min_read_quality", "NA")} ({depth:.0f}x)'),
                         ("Medaka model", pol.get("medaka_model", "NA")),
                         ("Length before -> after polishing (bp)",
                          f'{pol.get("draft_length", "NA")} -> {pol.get("polished_length", "NA")}')]
    for f in find("contigs_stats.tsv"):
        for line in open(f):
            c = line.rstrip("\n").split("\t")
            if len(c) >= 6 and not line.startswith(("#", "contig_id")) and c[0] != "final_mitogenome" and pol:
                mito_summary.append(("ONT source contig (pre-purge assembly), circularised by MitoHiFi", f"{c[0]}, {c[5]}"))
    if check_rows:
        bad = [r[0] for r in check_rows if r[4] != "OK"]
        mito_summary.append(("Protein-coding genes full length vs reference",
                             f"{len(check_rows) - len(bad)} of {len(check_rows)}"
                             + (f" (not OK: {', '.join(bad)}; probable sequence error)" if bad else "")))
    for f in find("mito_coverage.txt"):
        for line in open(f):
            c = line.split("\t")
            if not line.startswith("#") and len(c) >= 7:
                mito_summary.append((f"ONT reads on the mapping reference ({c[0]})",
                                     f"{c[3]} reads, {float(c[5]):.0f}% covered, mean depth {float(c[6]):.0f}x"))
    table("stage05_mito_summary", "Mitogenome: sequence and quality summary",
          "Source, length and protein-coding gene check of the mitogenome.",
          ["Item", "Value"], mito_summary, **PARENT_MITO)
    table("stage05_mito_gene_check", "Mitogenome: protein-coding gene check",
          "Each lifted CDS translated with the mitochondrial genetic code and compared with the reference "
          "protein. OK = same length and no internal stop.",
          ["Gene", "Protein length (aa)", "Reference length (aa)", "Internal stops", "Status"], check_rows,
          **PARENT_MITO)

lift = find("mitogenome.liftoff.gff3")
if lift:
    feats = gff_rows(lift[0])
    kind = {a.get("Parent"): t for t, s, e, st, a in feats if t in ("CDS", "tRNA", "rRNA")}
    prod = {a.get("Parent"): a.get("product", "") for t, s, e, st, a in feats if t in ("CDS", "tRNA", "rRNA")}
    parts = Counter(a.get("Parent") for t, s, e, st, a in feats if t in ("CDS", "exon"))
    rows = []
    for t, s, e, st, a in feats:
        if t != "gene":
            continue
        k = kind.get(a.get("ID"), "gene")
        rows.append((a.get("gene") or prod.get(a.get("ID")) or a.get("Name", "?"),
                     "protein-coding" if k == "CDS" else k, s, e, st, e - s + 1,
                     parts.get(a.get("ID"), 1), prod.get(a.get("ID"), ""), a.get("sequence_ID", "")))
    # distinct row names (two tRNA-Met)
    seen = Counter()
    named = []
    for r in sorted(rows, key=lambda r: r[2]):
        seen[r[0]] += 1
        named.append(((r[0] if seen[r[0]] == 1 else f"{r[0]} ({seen[r[0]]})"),) + r[1:])
    unmapped = [l.strip() for f in find("liftoff_unmapped_features.txt") for l in open(f) if l.strip()]
    table("stage05_mito_liftoff", "Mitogenome: annotation by liftover (Liftoff)",
          "Curated annotation of the reference mitogenome transferred onto the final sequence. Gene span includes introns. "
          + (f"Reference genes NOT transferred: {', '.join(unmapped)}." if unmapped else "All reference genes were transferred."),
          ["Gene", "Type", "Start", "End", "Strand", "Span (bp)", "Exons/CDS parts", "Product", "Identity to reference"],
          named, **PARENT_MITO)

mitos = find("mitogenome.mitos2.gff")
if mitos:
    feats = gff_rows(mitos[0])
    kind, parents, rows = {}, set(), []
    for t, s, e, st, a in feats:
        if t in ("mRNA", "tRNA", "rRNA"):
            kind[a.get("Parent")] = "protein-coding" if t == "mRNA" else t
            parents.add(a.get("ID"))
    for t, s, e, st, a in feats:
        if t in ("gene", "ncRNA_gene"):
            n = a.get("gene_id", a.get("ID", "?"))
            rows.append((n, kind.get(a.get("ID"), "gene"), s, e, st, e - s + 1,
                         "fragment or extra copy" if re.search(r"_\d+$", n) else ""))
        elif t == "exon" and a.get("Parent") not in parents:
            rows.append((a.get("Name", "?"), "fragment only", s, e, st, e - s + 1, "no gene model built by MITOS2"))
    seen = Counter()
    named = []
    for r in sorted(rows, key=lambda r: r[2]):
        seen[r[0]] += 1
        named.append(((r[0] if seen[r[0]] == 1 else f"{r[0]} ({seen[r[0]]})"),) + r[1:])
    table("stage05_mito_mitos2", "Mitogenome: de novo annotation (MITOS2)",
          "De novo annotation with MITOS2. 'fragment only' = hit without a gene model.",
          ["Feature", "Type", "Start", "End", "Strand", "Length (bp)", "Note"], named, **PARENT_MITO)

# ---- 06 repeats: RepeatMasker ---------------------------------------------------
PARENT_REPEATS = dict(
    parent_id="stage06_repeats", parent_name="06 Repeats",
    parent_description="Repeat and rRNA-gene annotation of the assembly.")
rows = []
for f in find("repeatmasker.tbl"):
    for line in open(f):
        m = re.match(r"\s*(bases masked):\s+(\d+) bp\s+\(\s*([\d.]+) %\)", line)
        if m:
            rows.append(("Total masked", "", m.group(2), m.group(3)))
            continue
        # top-level classes (no indent) and their main subclasses (3-space indent)
        m = re.match(r"^( {0}| {3})(\S[^:]*?):?\s+(\d+)\s+(\d+) bp\s+([\d.]+) %", line)
        if m:
            label = ("  " if m.group(1) else "") + m.group(2).strip()
            rows.append((label, m.group(3), m.group(4), m.group(5)))
        m2 = re.match(r"^(Total interspersed repeats):\s+(\d+) bp\s+([\d.]+) %", line)
        if m2:
            rows.append((m2.group(1), "", m2.group(2), m2.group(3)))
table("stage06_repeatmasker", "RepeatMasker",
      "Genome fraction per repeat class, masked with the combined RepeatModeler + EDTA library.",
      ["Repeat class", "Elements", "Length (bp)", "% of genome"], rows, **PARENT_REPEATS)

# ---- 06 rRNA: barrnap -----------------------------------------------------------
cnt, partial = Counter(), Counter()
for f in find("rrna.gff"):
    for line in open(f):
        if line.startswith("#") or "\t" not in line:
            continue
        attrs = line.rstrip("\n").split("\t")[8]
        name = re.search(r"Name=([^;]+)", attrs)
        n = name.group(1) if name else "rRNA"
        cnt[n] += 1
        if "partial" in attrs:
            partial[n] += 1
table("stage06_barrnap", "rRNA genes (barrnap)",
      "rRNA genes predicted by barrnap 0.9 (-k euk); not part of the gene annotation.",
      ["rRNA", "Predicted", "of which partial"], [(k, cnt[k], partial[k]) for k in sorted(cnt)],
      **PARENT_REPEATS)

# ---- 07 Structural annotation (BRAKER4 + supporting tools) ----------------------
# Everything in this stage - BRAKER4 itself plus the tools that feed or check it
# (STAR, StringTie, compleasm, OMArk, tRNAscan-SE, Infernal, GFF3 validation) -
# is grouped under one "07 Structural" nav entry so the report doesn't show a
# dozen flat top-level items for a single pipeline stage.
PARENT_STRUCTURAL = dict(
    parent_id="stage07_structural", parent_name="07 Structural",
    parent_description="BRAKER4 gene prediction and its supporting/QC tools.")
gss = find("gene_set_statistics.txt")
if gss:
    rows = []
    for line in open(gss[0]):
        m = re.match(r"^([A-Za-z][^:]+):\s+(.+)$", line.strip())
        if m:
            rows.append((m.group(1).strip(), m.group(2).strip()))
    fstats = [f for f in find("*.stats.json") if "annotate_results" in f]
    if fstats:
        pc = protein_coding(json.load(open(fstats[0])).get("annotation", {}))
        if pc:
            rows.append(("Final protein-coding genes (after merge/refinement)", pc))
    table("stage07_braker_geneset", "BRAKER4 raw gene set (before merge/refinement)", "BRAKER4 gene set (ETP mode) before merging and refinement. "
          "Last row: protein-coding genes in the final annotation.", ["Statistic", "Value"], rows,
          **PARENT_STRUCTURAL)
bs = find("busco_summary.txt")
busco_rows = {}
if bs:
    rows, section = [], None
    for line in open(bs[0]):
        s = line.strip()
        if s.startswith("---"):
            section = s.strip("- ")
        m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%", s)
        if m and section:
            rows.append((section, *m.groups()))
            busco_rows[section] = tuple(float(x) for x in m.groups())
    table("stage07_braker_busco", "BRAKER4: BUSCO (genome vs predicted proteome)",
          "BUSCO completeness of the assembly and of the BRAKER4 proteome.",
          ["Input", "Complete %", "Single %", "Duplicated %", "Fragmented %", "Missing %"], rows,
          **PARENT_STRUCTURAL)

# Completeness plot (BUSCO genome/proteome + compleasm genome/proteome): native
# 100%-stacked bargraph, same S/D/F/M categories/colours as the original PNG.
SDFM_CATS = [
    {"name": "Single", "color": "#2271B2"},
    {"name": "Duplicated", "color": "#9FC5E8"},
    {"name": "Fragmented", "color": "#F1C232"},
    {"name": "Missing", "color": "#CC0000"},
]
def compleasm_vals(path):
    if not os.path.exists(path):
        return None
    txt = open(path).read()
    v = {}
    for k, label in (("S", "Single"), ("D", "Duplicated"), ("F", "Fragmented"), ("M", "Missing")):
        m = re.search(rf"^{k}:([\d.]+)%", txt, re.M)
        v[label] = float(m.group(1)) if m else 0.0
    return v

completeness_data = {}
if "Genome Assembly" in busco_rows:
    c, s, d, f, m = busco_rows["Genome Assembly"]
    completeness_data["BUSCO genome"] = {"Single": s, "Duplicated": d, "Fragmented": f, "Missing": m}
if "Predicted Proteome" in busco_rows:
    c, s, d, f, m = busco_rows["Predicted Proteome"]
    completeness_data["BUSCO proteome"] = {"Single": s, "Duplicated": d, "Fragmented": f, "Missing": m}
cg = (find("compleasm_genome_out/summary.txt") or find("*compleasm_genome*/summary.txt"))
if cg:
    v = compleasm_vals(cg[0])
    if v:
        completeness_data["compleasm genome"] = v
cp = find("compleasm_summary.txt")
if cp:
    v = compleasm_vals(cp[0])
    if v:
        completeness_data["compleasm proteome (all isoforms)"] = v
plot("stage07_braker_completeness", "BRAKER4: completeness (BUSCO + compleasm)",
     "Genome and proteome completeness, BUSCO and compleasm (anthozoa_odb12). The compleasm proteome run uses all isoforms.",
     "bargraph", {"ylab": "% of BUSCOs", "stacking": "relative"}, completeness_data,
     cats=SDFM_CATS, **PARENT_STRUCTURAL)

# Transcript length distributions: computed directly from braker.gtf.gz (CDS
# length = sum of CDS feature lengths per transcript; genomic span = CDS-bounded
# range, matching how BRAKER's own gene_set_statistics.txt defines "genomic span").
# Rendered as a linegraph, not a bargraph: that's the same convention MultiQC's
# own FastQC module uses for its length-distribution histograms (see
# fastqc.seq_length_dist_plot) - a fine-grained line reads as a histogram
# silhouette, whereas bargraph draws gapped categorical bars that don't.
def histogram_points(values, bin_width, cap):
    """[(bin_start, count), ...], numeric x, sorted. Values above `cap` are
    pooled into the last point (so the x-axis doesn't stretch to one outlier).
    Returned as pairs, not a {x: y} dict - see the note by cum_data above."""
    bins = Counter()
    for v in values:
        bins[min((v // bin_width) * bin_width, cap)] += 1
    return [(x, bins[x]) for x in sorted(bins)]

gtf = find("braker.gtf.gz") or find("*.gtf.gz")
if gtf:
    cds_len, cds_bounds = Counter(), {}
    with gzip.open(gtf[0], "rt") as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "CDS":
                continue
            m = re.search(r'transcript_id "([^"]+)"', c[8])
            if not m:
                continue
            tid = m.group(1)
            s, e = int(c[3]), int(c[4])
            cds_len[tid] += e - s + 1
            if tid not in cds_bounds:
                cds_bounds[tid] = [s, e]
            else:
                cds_bounds[tid][0] = min(cds_bounds[tid][0], s)
                cds_bounds[tid][1] = max(cds_bounds[tid][1], e)
    span_len = {t: (e - s + 1) for t, (s, e) in cds_bounds.items()}
    cds_hist = histogram_points(cds_len.values(), bin_width=200, cap=10000)
    span_hist = histogram_points(span_len.values(), bin_width=1000, cap=50000)
    plot("stage07_braker_transcript_lengths_cds", "BRAKER4 raw gene set: transcript lengths (CDS)",
         "CDS length per transcript of the raw BRAKER4 set, 200 bp bins (&ge;10,000 bp pooled into the last point).",
         "linegraph", {"xlab": "CDS length (bp)", "ylab": "Transcripts", "ymin": 0,
                        "showlegend": False},
         {"Transcripts": cds_hist}, **PARENT_STRUCTURAL)
    plot("stage07_braker_transcript_lengths_span", "BRAKER4 raw gene set: transcript lengths (genomic span)",
         "Genomic span per transcript of the raw BRAKER4 set, 1,000 bp bins (&ge;50,000 bp pooled into the last point).",
         "linegraph", {"xlab": "Genomic span (bp)", "ylab": "Transcripts", "ymin": 0,
                        "showlegend": False},
         {"Transcripts": span_hist}, **PARENT_STRUCTURAL)

# ---- 07 STAR alignment QC flag ---------------------------------------------------
# The STAR module's own general-stats numbers have no pass/fail signal attached,
# so a bad RNA-seq alignment rate (which directly weakens BRAKER's hints and the
# StringTie assembly above) is easy to miss. Flag it explicitly per sample.
star_rows = []
for f in find("*_Log.final.out"):
    kv = {}
    for line in open(f):
        if "|" in line:
            k, v = line.split("|", 1)
            kv[k.strip()] = v.strip()
    uniq = kv.get("Uniquely mapped reads %", "").rstrip("%")
    multi = kv.get("% of reads mapped to multiple loci", "").rstrip("%")
    try:
        uniq_f = float(uniq)
    except ValueError:
        continue
    status = "FAIL" if uniq_f < 50 else "WARN" if uniq_f < 80 else "PASS"
    name = os.path.basename(f)
    m = re.search(r"Por_([A-Za-z0-9]+_\d+_\d+)", name)
    label = f"RNA-seq {m.group(1)}" if m else stem(f, "_Log.final.out")
    star_rows.append((label, uniq, multi, status))
if star_rows:
    star_rows.sort(key=lambda r: float(r[1]))
    table("stage07_star_qc_flag", "STAR: RNA-seq alignment QC flag",
          "Uniquely-mapped read rate per RNA-seq sample (FAIL &lt;50%, WARN &lt;80%, PASS &ge;80%).",
          ["Sample", "% Uniquely mapped", "% Multi-mapped", "Status"], star_rows,
          **PARENT_STRUCTURAL)

# ---- 07 tRNAscan-SE -------------------------------------------------------------
def trnascan_counts(path):
    trna, pseudo = Counter(), Counter()
    for line in open(path):
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 10 or not cols[1].strip().isdigit():
            continue
        iso = cols[4].strip()
        (pseudo if "pseudo" in line else trna)[iso] += 1
    return trna, pseudo

trna, pseudo = Counter(), Counter()
for f in find("trnascan.out"):
    t, p = trnascan_counts(f)
    trna.update(t)
    pseudo.update(p)
hc = None
for f in find("trna_hc.out"):
    hc, _ = trnascan_counts(f)
fin_trna = None
fg = final_gff()
if fg and hc is not None:   # refined runs only; older funannotate GFFs have product=None on tRNAs
    fin_trna = Counter()
    for line in open(fg):
        c = line.rstrip("\n").split("\t")
        if len(c) > 8 and c[2] == "tRNA":
            m = re.search(r"isotype=([^;]+)", c[8]) or re.search(r"product=tRNA-([^;]+)", c[8])
            fin_trna[m.group(1) if m else "?"] += 1
if trna or pseudo:
    keys = sorted(set(trna) | set(pseudo) | set(hc or {}) | set(fin_trna or {}))
    header = ["Isotype", "Raw tRNAscan-SE tRNAs", "Raw pseudogenes"]
    if hc is not None:
        header.append("High-confidence")
    if fin_trna is not None:
        header.append("In final annotation")
    rows = []
    for k in keys + ["Total"]:
        get = (lambda c: sum(c.values())) if k == "Total" else (lambda c: c[k])
        r = [k, get(trna), get(pseudo)]
        if hc is not None:
            r.append(get(hc))
        if fin_trna is not None:
            r.append(get(fin_trna))
        rows.append(r)
    if fin_trna is not None:
        desc = ("tRNAs per isotype: tRNAscan-SE calls and pseudogenes, the EukHighConfidenceFilter set, and "
                "the tRNAs in the final annotation.")
    else:
        desc = "tRNAscan-SE calls and pseudogenes per isotype."
    table("stage07_trnascan", "tRNA genes (tRNAscan-SE)", desc, header, rows, **PARENT_STRUCTURAL)

# ---- 07 ncRNA: Infernal / Rfam --------------------------------------------------
fam, typ = Counter(), Counter()
for f in find("cmscan.tblout"):
    for line in open(f):
        if line.startswith("#"):
            continue
        c = line.split()
        # --fmt 2 (1-based): col 19 is the inclusion flag ('!'), col 20 is 'olp'
        # ('*' no overlap, '^' best of an overlap group, '=' overlapped by a better hit).
        if len(c) < 21 or c[18] != "!" or c[19] not in ("*", "^"):
            continue
        fam[c[1]] += 1
for k, v in fam.items():
    kl = k.lower()
    t = ("tRNA" if kl.startswith(("trna", "selenocysteine_trna")) else
         "rRNA" if "rrna" in kl or kl.startswith(("ssu", "lsu", "5s", "5_8s")) else
         "snoRNA" if kl.startswith(("snor", "snord", "snora", "scarna", "u3", "u8")) else
         "snRNA/spliceosomal" if re.match(r"^u\d|u1|u2|u4|u5|u6|u11|u12", kl) else
         "miRNA" if kl.startswith("mir") else
         "other")
    typ[t] += v
table("stage07_infernal_types", "Rfam hits (raw cmscan), by type",
      "All significant, non-overlapping Rfam hits (cmscan --cut_ga); type from the Rfam family name.",
      ["Type", "Hits"], sorted(typ.items(), key=lambda x: -x[1]), **PARENT_STRUCTURAL)
table("stage07_infernal_families", "Rfam hits (raw cmscan), top families",
      "The 30 most frequent Rfam families among all significant hits.",
      ["Rfam family", "Hits"], fam.most_common(30), **PARENT_STRUCTURAL)
rg = refined_gff()
if rg:
    cls = Counter()
    for line in open(rg):
        c = line.rstrip("\n").split("\t")
        if len(c) > 8 and c[2] in ("ncRNA", "rRNA") and "RFAM:" in c[8]:
            if c[2] == "rRNA":
                m = re.search(r"Name=([^;]+)", c[8])
                cls["rRNA: " + (m.group(1) if m else "?")] += 1
            else:
                m = re.search(r"ncRNA_class=([^;]+)", c[8])
                cls[m.group(1) if m else "other"] += 1
    table("stage07_ncrna_final", "ncRNAs in the final annotation, by class",
          "Rfam-derived ncRNA and rRNA genes in the final annotation.",
          ["Class", "Genes"], sorted(cls.items(), key=lambda x: -x[1]), **PARENT_STRUCTURAL)

# ---- 07 StringTie, compleasm, OMArk, GFF3 validation -----------------------------
# When a refined GFF exists but was made after the pipeline (the v1.0 release workflow), OMArk and
# StringTie still describe the pre-refinement model set; the refinement does not change any CDS, so
# the proteins are identical. The pipeline refines BEFORE those steps and sets MQC_REFINED_UPSTREAM=1.
PRE_REFINE_NOTE = (" Gene models before refinement."
                   if refined_gff() and os.environ.get("MQC_REFINED_UPSTREAM") != "1" else "")
for f in find("stringtie_stats.txt"):
    rows = [tuple(l.split("\t")[:2]) for l in open(f).read().splitlines() if l and not l.startswith("#") and "\t" in l]
    table("stage07_stringtie", "StringTie transcript assembly",
          "Reference-guided transcript assembly from the RNA-seq alignments." + PRE_REFINE_NOTE, ["Metric", "Value"],
          rows, **PARENT_STRUCTURAL)
for f in find("summary.txt"):
    if "compleasm" not in f:
        continue
    txt = open(f).read()
    rows = [(m.group(1), m.group(2), m.group(3)) for m in re.finditer(r"^([SDFIM]):([\d.]+)%, (\d+)", txt, re.M)]
    lin = re.search(r"lineage: (\S+)", txt)
    table("stage07_compleasm", "compleasm (genome completeness)",
          f"compleasm on the masked assembly ({lin.group(1) if lin else ''}). S single, D duplicated, F fragmented, I incomplete, M missing.",
          ["Category", "%", "Count"], rows, **PARENT_STRUCTURAL)
    break
for f in find("*.sum"):
    if "omark" not in f:
        continue
    txt = open(f).read()
    rows = []
    m = re.search(r"^S:([\d.]+)%,D:([\d.]+)%\[U:([\d.]+)%,E:([\d.]+)%\],M:([\d.]+)%", txt, re.M)
    if m:
        rows += [("Single (conserved HOGs)", m.group(1)), ("Duplicated", m.group(2)),
                 ("  of which unexpected", m.group(3)), ("  of which expected", m.group(4)),
                 ("Missing", m.group(5))]
    m = re.search(r"^A:([\d.]+)%\[P:([\d.]+)%,F:([\d.]+)%\],I:([\d.]+)%.*?C:([\d.]+)%.*?U:([\d.]+)%", txt, re.M)
    if m:
        rows += [("Taxonomically consistent", m.group(1)), ("Inconsistent", m.group(4)),
                 ("Likely contamination", m.group(5)), ("Unknown", m.group(6))]
    n = re.search(r"there are (\d+) proteins", txt)
    clade = re.search(r"selected clade was (.+)", txt)
    table("stage07_omark", "OMArk (proteome completeness and consistency)",
          f"{n.group(1) if n else ''} proteins (one isoform per gene), clade {clade.group(1).strip() if clade else ''}.{PRE_REFINE_NOTE}",
          ["Metric", "% "], rows, **PARENT_STRUCTURAL)
    break
def gt_table(f, sec_id, name, description):
    # Bucket gt gff3validator messages by kind (full text, IDs/numbers masked),
    # counting how many messages and how many distinct seqids each kind covers.
    kinds, seqids, example = Counter(), defaultdict(set), {}
    for l in open(f).read().splitlines():
        if not (l.startswith("warning") or l.startswith("error")):
            continue
        level, msg = l.split(":", 1)
        msg = msg.strip()
        sid = re.search(r'seqid "([^"]+)"', msg)
        kind = level.upper() + ": " + re.sub(r'"[^"]*"', '"…"', re.sub(r"line \d+", "line …", msg))
        kind = re.sub(r'in file "…"', "", kind).strip()
        kinds[kind] += 1
        if sid:
            seqids[kind].add(sid.group(1))
        example.setdefault(kind, msg)
    SHORT = {
        'has not been previously introduced with a "…" line':
            "Missing ##sequence-region pragma",
    }
    rows = []
    for kind, n in kinds.most_common():
        label = next((v for k, v in SHORT.items() if k in kind), kind)
        rows.append((label, n, len(seqids[kind]) or "", example[kind]))
    if not rows:
        rows = [("No warnings or errors", 0, "", "")]
    table(sec_id, name, description, ["Message kind", "Count", "Seqids affected", "Example"], rows,
          **PARENT_STRUCTURAL)

refined_val = find("gt_refined.txt")
pre_val = find("gff3_validation.txt")
if refined_val:
    gt_table(refined_val[0], "stage07_gff3_validation", "GFF3 validation (GenomeTools)",
             "gt gff3validator messages on the released annotation, by kind, with the number of seqids affected and one example.")
    if pre_val:
        gt_table(pre_val[0], "stage07_gff3_validation_prerefine", "GFF3 validation, pre-refinement (GenomeTools)",
                 "gt gff3validator messages on the merged annotation before refinement.")
elif pre_val:
    gt_table(pre_val[0], "stage07_gff3_validation", "GFF3 validation (GenomeTools)",
             "gt gff3validator messages on the merged annotation, bucketed by kind, with the number of seqids "
             "affected and one full example message.")

# ---- 08 functional annotation: funannotate --------------------------------------
PARENT_FUNCTIONAL = dict(
    parent_id="stage08_functional", parent_name="08 Functional",
    parent_description="Functional annotation (funannotate) and membrane/signal peptide prediction (Phobius).")
stats = [f for f in find("*.stats.json") if "annotate_results" in f] or find("annotation_stats.json")
if stats:
    a = json.load(open(stats[0])).get("annotation", {})
    tl = a.get("transcript-level", {})
    fn = tl.get("functional", {})
    rows = [("Protein-coding genes", protein_coding(a)), ("All genes (incl. tRNA/ncRNA/rRNA)", a.get("genes")),
            ("mRNAs", a.get("mRNA")), ("tRNAs", a.get("tRNA")),
            ("ncRNAs", a.get("ncRNA")), ("rRNAs", a.get("rRNA")),
            ("Genes with a gene symbol", a.get("common_name")),
            ("Mean protein length (aa)", tl.get("avg_protein_length")),
            ("Complete CDS", tl.get("CDS_complete")), ("Single-exon transcripts", tl.get("single_exon_transcript"))]
    rows += [(f"Functional: {k}", v) for k, v in fn.items()]
    table("stage08_funannotate", "Annotation summary (funannotate)",
          "funannotate annotate: gene counts and functional annotation (eggNOG, InterProScan, Phobius, Pfam, dbCAN, MEROPS, SwissProt, BUSCO).",
          ["Metric", "Value"], [r for r in rows if r[1] is not None], **PARENT_FUNCTIONAL)

# ---- 08 Phobius -----------------------------------------------------------------
for f in find("phobius_results.txt"):
    tm = sp = n = 0
    for line in open(f):
        c = line.split()
        if len(c) < 3 or c[0].startswith("SEQ") or not c[1].isdigit():
            continue
        n += 1
        tm += int(c[1]) > 0
        sp += c[2] == "Y"
    table("stage08_phobius", "Transmembrane and signal peptides (Phobius)",
          "Proteins with at least one predicted transmembrane helix or a signal peptide.",
          ["Class", "Proteins"], [("Proteins scanned", n), ("With TM helix", tm), ("With signal peptide", sp)],
          **PARENT_FUNCTIONAL)
    break

# ---- 10 NCBI validation (table2asn .val) ----------------------------------------
for f in find("*.val")[:1]:
    sev = Counter()
    for line in open(f, errors="replace"):
        m = re.match(r"(Info|Warning|Error|Reject|Fatal|Critical)\S*: valid \[([A-Za-z0-9_.]+)\]", line, re.I)
        if m:
            sev[(m.group(1).capitalize(), m.group(2))] += 1
    rows = [(f"{s_} {code}", n) for (s_, code), n in sorted(sev.items(), key=lambda kv: (-kv[1], kv[0]))]
    table("stage10_table2asn_val", "NCBI validation (table2asn)",
          "Validator messages in the submission file, by severity and code. The QC gate fails on any "
          "Error, Reject, Fatal or Critical.", ["Severity and code", "Count"], rows)

# ---- 00 NanoPlot stats: give raw/trimmed files distinct names -------------------
# Both NanoPlot runs write NanoStats.txt, so MultiQC saw one sample name twice and
# silently kept only one ("Duplicate sample data found! Overwriting: NanoStats").
# Copy each under a name carrying its read set; the originals are ignored via
# fn_ignore_paths in multiqc_config.yaml.
for f in find("NanoStats*.txt"):
    parts = f.split(os.sep)
    src = next((p for p in parts if p.startswith("qc_nanoplot_")), None)
    if src:
        tag = src.replace("qc_nanoplot_", "ONT_")
        shutil.copyfile(f, os.path.join(OUT, f"{tag}.{os.path.basename(f)}"))

# ---- Pipeline summary (top of report) -------------------------------------------
# Aggregates values already computed elsewhere; nothing is recomputed.
def gate_info():
    for f in find("polish_qc_gate_mqc.tsv"):
        txt = open(f).read()
        fin = re.search(r"Final assembly: (\w+)", txt)
        verdicts = {}
        for line in txt.splitlines():
            c = line.split("\t")
            if len(c) == 5 and c[0] != "metric":
                verdicts[c[0]] = (c[1], c[2], c[4])
        return (fin.group(1).lower() if fin else None), verdicts
    for f in find("polish_qc_report.tsv"):   # older runs: no mqc table
        verdicts = {}
        for line in open(f).read().splitlines()[1:]:
            c = line.split("\t")
            if len(c) == 5:
                verdicts[c[0]] = (c[1], c[2], c[4])
        fail = any(v[2] == "FAIL" for v in verdicts.values())
        return ("unpolished" if fail else "polished"), verdicts
    return None, {}

def busco_c(stage):
    for f in find(f"{stage}.*short_summary*.txt"):
        m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%", open(f).read())
        if m:
            return m.groups()
    return None

def gfa(stage):
    for f in find(f"{stage}.gfastats.txt"):
        kv = {}
        for line in open(f):
            if ":" in line:
                k, v = line.split(":", 1)
                kv[k.strip()] = v.strip()
        return kv
    return {}

decision, verdicts = gate_info()
if decision:
    final_busco_stage = "final_unpolished" if decision == "unpolished" else "postpolish"
    # gfastats is run on the gate inputs: the unpolished final product is the
    # prepolish (filtered) assembly after funannotate clean/sort.
    final_gfa_stage = "prepolish" if decision == "unpolished" else "postpolish"
    rows = [("Final assembly carried forward",
             "unpolished" if decision == "unpolished" else "polished (Medaka)", "")]
    g = gfa(final_gfa_stage)
    note = " (before funannotate clean/sort)" if decision == "unpolished" else ""
    if g:
        rows += [("Genome size (bp)" + note, g.get("Total scaffold length", "NA"), ""),
                 ("Scaffolds" + note, g.get("# scaffolds", "NA"), ""),
                 ("Scaffold N50 (bp)", g.get("Scaffold N50", "NA"), ""),
                 ("Contig N50 (bp)", g.get("Contig N50", "NA"), ""),
                 ("GC %", g.get("GC content %", "NA"), "")]
    b = busco_c(final_busco_stage)
    if b:
        rows.append(("BUSCO complete % (final assembly)", b[0], ""))
        rows.append(("BUSCO missing % (final assembly)", b[4], ""))
    for metric, label in (("merqury_qv", "Polish QC gate: Merqury QV (pre -> post)"),
                          ("busco_completeness", "Polish QC gate: BUSCO C% (pre -> post)")):
        if metric in verdicts:
            pre, post, v = verdicts[metric]
            rows.append((label, f"{pre} -> {post}", v))
    rows.append(("Polish QC gate: decision",
                 "unpolished" if decision == "unpolished" else "polished",
                 "FAIL" if decision == "unpolished" else "PASS"))
    stats = [f for f in find("*.stats.json") if "annotate_results" in f]
    if stats:
        a = json.load(open(stats[0])).get("annotation", {})
        rows += [("Protein-coding genes (final annotation)", protein_coding(a) or "NA", ""),
                 ("All genes incl. tRNA/ncRNA/rRNA (final annotation)", a.get("genes", "NA"), ""),
                 ("mRNAs / tRNAs / ncRNAs / rRNAs",
                  f'{a.get("mRNA", "NA")} / {a.get("tRNA", "NA")} / {a.get("ncRNA", "NA")} / {a.get("rRNA", "NA")}', "")]
    d = dict(mito_summary)
    if d:
        rows.append(("Mitogenome (separate from the nuclear assembly)",
                     f'{d.get("Length (bp)", "NA")} bp; {d.get("Sequence source", "NA")}; protein-coding genes '
                     f'full length: {d.get("Protein-coding genes full length vs reference", "not checked")}',
                     "" if "Protein-coding genes full length vs reference" not in d
                     else ("FAIL" if "not OK" in d["Protein-coding genes full length vs reference"] else "PASS")))
    table("pipeline_summary", "Pipeline summary",
          "Final assembly and annotation, with QC gate results.",
          ["Item", "Value", "Gate status"], rows)

print("\n".join(sorted(os.listdir(OUT))))
