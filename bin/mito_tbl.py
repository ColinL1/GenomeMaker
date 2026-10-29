#!/usr/bin/env python3
"""Mitogenome for the genome submission: FASTA with NCBI source modifiers and a 5-column feature table
from the Liftoff annotation (curated RefSeq reference), for table2asn together with the nuclear scaffolds.

NCBI: organelle sequences go in the same submission as the genome, labelled in the definition line
([location=mitochondrion]); the mitochondrial genetic code comes from the organism's taxonomy.

usage: mito_tbl.py mitogenome.fasta mitogenome.liftoff.gff3 mito_gene_check.tsv LOCUS_PREFIX out.fsa out.tbl
The sequence ID is <LOCUS_PREFIX>_mtDNA (not 'MT' or 'mitogenome': NCBI's discrepancy report flags both
as SUSPICIOUS_SEQUENCE_ID). Genes get locus tags <PREFIX>_MT01, _MT02, ... in genome order. CDS split by
introns become joins.

Source modifiers: [topology=circular] because both mitogenome routes (stage 05) only pass on a circular
molecule (GetOrganelle 'circular genome', MitoHiFi circularised), and MITO_ROTATE keeps that single
record. [completeness=complete] only when every reference protein is found full length
(mito_gene_check.tsv, all 'OK'); otherwise the modifier is left out and the reason printed.
"""
import collections
import os
import re
import sys

fa, gff, check, prefix, out_fa, out_tbl = sys.argv[1:7]
sid = f'{prefix}_mtDNA'

recs = [l for l in open(fa) if l.startswith('>')]
if len(recs) != 1:
    sys.exit(f'ERROR: {fa} has {len(recs)} records; the submission expects one circular mitogenome')
seq = ''.join(l.strip() for l in open(fa) if not l.startswith('>'))
st = [l.rstrip('\n').split('\t') for l in open(check)][1:] if not os.path.basename(check).startswith('NO_') else []
not_ok = [f'{r[0]} {r[-1]}' for r in st if r[-1] != 'OK']
complete = bool(st) and not not_ok


def attrs(s):
    return dict(kv.split('=', 1) for kv in s.strip().rstrip(';').split(';') if '=' in kv)


genes, kids = {}, collections.defaultdict(list)
for line in open(gff):
    if line.startswith('#') or not line.strip():
        continue
    f = line.rstrip('\n').split('\t')
    at = attrs(f[8])
    if f[2] == 'gene':
        genes[at['ID']] = dict(start=int(f[3]), end=int(f[4]), strand=f[6], name=at.get('gene', ''),
                               biotype=at.get('gene_biotype', ''))
    elif f[2] in ('CDS', 'rRNA', 'tRNA') and 'Parent' in at:
        kids[at['Parent']].append((f[2], int(f[3]), int(f[4]), f[6], at.get('product', '')))


def loc(a, b, strand):
    return (a, b) if strand == '+' else (b, a)


lines = [f'>Feature {sid}']
for i, (gid, g) in enumerate(sorted(genes.items(), key=lambda x: x[1]['start']), 1):
    a, b = loc(g['start'], g['end'], g['strand'])
    lines.append(f'{a}\t{b}\tgene')
    if g['biotype'] == 'protein_coding' and g['name']:
        lines.append(f'\t\t\tgene\t{g["name"]}')
    lines.append(f'\t\t\tlocus_tag\t{prefix}_MT{i:02d}')
    feats = kids[gid]
    types = {t for t, *_ in feats}
    for t in ('CDS', 'rRNA', 'tRNA'):
        segs = sorted((x, y, s, p) for tt, x, y, s, p in feats if tt == t)
        if not segs:
            continue
        if segs[0][2] == '-':
            segs = segs[::-1]
        for k, (x, y, s, p) in enumerate(segs):
            c = loc(x, y, s)
            lines.append(f'{c[0]}\t{c[1]}' + (f'\t{t}' if k == 0 else ''))
        lines.append(f'\t\t\tproduct\t{segs[0][3]}')
        if t == 'CDS':   # NCBI requires a protein_id on every CDS (same gnl|ncbi| style as funannotate's)
            lines.append(f'\t\t\tprotein_id\tgnl|ncbi|{prefix}_MT{i:02d}')
    if not types:
        sys.exit(f'ERROR: gene {gid} has no CDS/rRNA/tRNA child')
with open(out_tbl, 'w') as o:
    o.write('\n'.join(lines) + '\n')
with open(out_fa, 'w') as o:
    o.write(f'>{sid} [location=mitochondrion]' + (' [completeness=complete]' if complete else '') + ' [topology=circular]\n')
    for i in range(0, len(seq), 80):
        o.write(seq[i:i + 80] + '\n')
print(f'{len(genes)} genes, {len(seq)} bp; completeness ' +
      ('complete' if complete else 'not stated (' + ('; '.join(not_ok) if not_ok else 'no gene check') + ')'))
