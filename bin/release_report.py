#!/usr/bin/env python3
"""README.md and summary.tsv for the release, one protein-coding number everywhere.

The protein-coding count is genes - tRNA - rRNA - ncRNA from funannotate's stats.json and is
checked against the number of genes with an mRNA in the release GFF3; the script fails if they differ.

usage: release_report.py --stats F.stats.json --release-gff R.gff3
         --summary summary.tsv --species 'Genus sp.' --prefix X --outdir DIR
         [--qc annotation_qc.tsv] [--version 1.0] [--isolate X] [--title T]
         [--methods methods.md] [--notes known_limitations.md]
The numbers in the README come from summary.tsv (release_summary.py). Optional files may be given as NO_* placeholders.
"""
import argparse
import json
import re
import sys

ap = argparse.ArgumentParser()
for k in ('stats', 'release-gff', 'qc', 'species', 'prefix', 'outdir', 'methods', 'notes'):
    ap.add_argument('--' + k)
ap.add_argument('--title', default='', help='README title; default "<species> genome <version>"')
ap.add_argument('--version', default='1.0')
ap.add_argument('--isolate', default='')
ap.add_argument('--summary', required=True, help='summary.tsv from release_summary.py (section, metric, value)')
a = ap.parse_args()


def usable(p):
    return p and not p.split('/')[-1].startswith('NO_')


an = json.load(open(a.stats))['annotation']
pc = an['genes'] - an.get('tRNA', 0) - an.get('rRNA', 0) - an.get('ncRNA', 0)
gene_types = {}
for line in open(a.release_gff):
    f = line.rstrip('\n').split('\t')
    if len(f) > 8 and f[2] in ('mRNA', 'tRNA', 'ncRNA', 'rRNA'):
        p = re.search(r'Parent=([^;]+)', f[8]).group(1)
        gene_types[p] = f[2]
pc_gff = sum(1 for t in gene_types.values() if t == 'mRNA')
if pc != pc_gff:
    sys.exit(f'ERROR: protein-coding genes differ: stats.json {pc}, release GFF3 {pc_gff}')
qc = [l.rstrip('\n').split('\t') for l in open(a.qc)][1:] if usable(a.qc) else []
qc_fail = [r for r in qc if r[2] == 'FAIL']
qc_warn = [r for r in qc if r[2] == 'WARN']

uni = [l.rstrip('\n').split('\t') for l in open(a.summary)][1:]
pc_sum = next(int(str(v).replace(',', '')) for sec, m, v in uni if m == 'Protein-coding genes')
if pc_sum != pc:
    sys.exit(f'ERROR: protein-coding genes differ: summary.tsv {pc_sum}, stats.json {pc}')
counts_md = '\n'.join(f'| {sec} | {m} | {v} |' for sec, m, v in uni)

P = a.prefix
title = a.title or f"{a.species} genome {a.version}" + (f" (isolate {a.isolate})" if a.isolate else '')
methods = open(a.methods).read().strip() if usable(a.methods) else ''
notes = open(a.notes).read().strip() if usable(a.notes) else ''
readme = f"""# {title}

## Summary

| Section | Metric | Value |
|---|---|---|
""" + counts_md + f"""

Protein-coding genes: **{pc:,}** (genes minus tRNA, rRNA and ncRNA genes in funannotate's stats.json;
the same number as the genes with an mRNA in `annotation/{P}.gff3`).

## Contents

| Folder | Files |
|---|---|
| `genome/` | `{P}.scaffolds.fa` (soft-masked, screened), `scaffold_name_map.tsv` |
| `annotation/` | `{P}.gff3` (canonical; with evidence/haplotig/gap tags), proteins, CDS and mRNA transcripts, `annotations.txt`, `.gbk` (GenBank flatfile of the `.sqn`, mitogenome included); GFF3, `annotations.txt` and `.gbk` carry the same products, EC numbers and gene spans as the submission |
| `annotation/extra/` | `{P}.hc.gff3`, `{P}.hc.proteins.fa`, `trna_pseudogenes.gff3`, `splits.tsv`, gap-fusion, gap-misorder and haplotig candidate lists, `id_map.tsv`, `rrna_decisions.tsv` (rRNA fragments: partial ends at a sequence end or gap, other fragments removed) |
| `repeats/` | `repeats.gff3`, `repeatmasker.tbl`, `repeatmasker.classified.out`, `repeats_classified.lib` |
| `ncbi/` | table2asn `.sqn`, `.val`, `.dr` (the mitogenome, when present, is in the same `.sqn` as [location=mitochondrion]) |
| `mitogenome/` | mitogenome sequence and annotation |
| `qc/` | `annotation_qc.tsv`, `refine_report.tsv`, screening reports, BUSCO/compleasm, symbiont BUSCO, Merqury, `multiqc_report.html` |
""" + (f"""
## Methods

{methods}
""" if methods else '') + (f"""
## Known limitations

{notes}
""" if notes else '')
if qc_warn or qc_fail:
    readme += '\n## QC gate notes\n\n' + '\n'.join(f'- {r[0]}: {r[2]} ({r[3]}; {r[4]})' for r in qc_fail + qc_warn) + '\n'
open(f'{a.outdir}/README.md', 'w').write(readme)
print(f'protein-coding genes {pc}; summary rows {len(uni)}')
