#!/usr/bin/env python3
"""Remove the contigs MitoHiFi identified as mitochondrial from the nuclear assembly.

MitoHiFi (contigs mode, stage 05) lists the contigs that match the reference mitogenome
(contigs_filtering/contigs_ids.txt). A listed contig longer than MAX_RATIO x the reference
mitogenome is kept: it is a nuclear contig carrying a mitochondrial insertion, not the mitogenome.

usage: remove_mito_contigs.py assembly.fa contigs_ids.txt contigs_stats.tsv mito_ref.fa MAX_RATIO out.fa report.tsv
       (contigs_ids.txt / contigs_stats.tsv may be NO_* placeholders: nothing is removed)
report.tsv: contig, length, circular, genes, action (removed | kept_too_long | not_in_assembly)
"""
import os
import sys

asm_p, ids_p, stats_p, ref_p, ratio, out_p, rep_p = sys.argv[1:8]
ratio = float(ratio)


def usable(p):
    return os.path.exists(p) and not os.path.basename(p).startswith('NO_') and os.path.getsize(p) > 0


def records(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith('>'):
            if name is not None:
                yield name, ''.join(buf)
            name, buf = line[1:].rstrip('\n'), []
        else:
            buf.append(line.strip())
    if name is not None:
        yield name, ''.join(buf)


ids = []
if usable(ids_p):
    ids = [l.split()[0] for l in open(ids_p) if l.strip()]
stats = {}
if usable(stats_p):
    for l in open(stats_p):
        f = l.rstrip('\n').split('\t')
        if len(f) >= 6 and not l.startswith(('#', 'contig_id')):
            stats[f[0]] = (f[5], f[4])   # was_circular, number_of_genes
ref_len = sum(len(s) for _, s in records(ref_p))
max_len = ratio * ref_len

lens = {}
removed = set()
for header, seq in records(asm_p):
    lens[header.split()[0]] = len(seq)
for c in ids:
    if c in lens and lens[c] <= max_len:
        removed.add(c)

with open(out_p, 'w') as o:
    for header, seq in records(asm_p):
        if header.split()[0] in removed:
            continue
        o.write(f'>{header}\n')
        for i in range(0, len(seq), 80):
            o.write(seq[i:i + 80] + '\n')

with open(rep_p, 'w') as r:
    r.write('contig\tlength\tcircular\tgenes\taction\n')
    for c in ids:
        circ, genes = stats.get(c, ('NA', 'NA'))
        if c not in lens:
            action = 'not_in_assembly'
        elif c in removed:
            action = 'removed'
        else:
            action = f'kept_too_long (> {ratio:g} x {ref_len} bp reference)'
        r.write(f'{c}\t{lens.get(c, "NA")}\t{circ}\t{genes}\t{action}\n')

print(f'MitoHiFi listed {len(ids)} contig(s); removed {len(removed)} '
      f'({sum(lens[c] for c in removed)} bp) from the nuclear assembly')
