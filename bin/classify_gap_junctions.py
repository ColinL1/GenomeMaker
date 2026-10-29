#!/usr/bin/env python3
"""WP6 classifier: label each gap junction from DIAMOND hits of its two flanking segments.

Classes (as in the 2026-10-16 review): same_protein_collinear, same_protein_noncollinear,
different_proteins, one_side_only, no_hit, segment_lt_30aa.

Tags:
  possible_gap_fusion   class different_proteins at an unsupported junction
  possible_gap_misorder class same_protein_noncollinear at an unsupported junction

Strict split rule (Plan 1, WP6). A junction is split only if ALL hold:
  1. support == unsupported (no RNA or protein intron hint);
  2. in at least MIN_SPECIES proteomes the two halves' best hits are different proteins;
  3. in those proteomes each half covers >= MIN_COV of its best reference protein;
  4. in those proteomes neither best reference protein is hit by both halves.
Conditions 2-4 are evaluated together per proteome; a proteome counts only if all three hold.

Reference protein IDs must be prefixed '<species>|' (one prefix per proteome).

usage: classify_gap_junctions.py junctions.tsv hits.tsv out_tags.tsv [min_species=2] [min_cov=0.6]
hits.tsv: DIAMOND -f 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen
"""
import collections
import re
import sys

junc_path, hits_path, out_path = sys.argv[1:4]
MIN_SPECIES = int(sys.argv[4]) if len(sys.argv) > 4 else 2
MIN_COV = float(sys.argv[5]) if len(sys.argv) > 5 else 0.6

# segment -> subject gene -> (sstart, send, bits, slen); best-bitscore HSP per subject.
# Subjects are compared at the gene level: isoforms '<gene>-T<n>' of one reference gene
# are the same protein for this test.
hits = collections.defaultdict(dict)
species = set()
for line in open(hits_path):
    f = line.rstrip('\n').split('\t')
    q, s, bits = f[0], re.sub(r'-T\d+$', '', f[1]), float(f[9])
    species.add(s.split('|', 1)[0])
    if s not in hits[q] or bits > hits[q][s][2]:
        hits[q][s] = (int(f[6]), int(f[7]), bits, int(f[11]))


def classify(left, right):
    L, R = hits.get(left), hits.get(right)
    if not L and not R:
        return 'no_hit'
    if not L or not R:
        return 'one_side_only'
    shared = L.keys() & R.keys()
    if not shared:
        return 'different_proteins'
    if any(L[s][0] < R[s][0] + 10 for s in shared):
        return 'same_protein_collinear'
    return 'same_protein_noncollinear'


def split_species(left, right):
    """Proteomes in which conditions 2-4 all hold; with the best hits used."""
    L, R = hits.get(left, {}), hits.get(right, {})
    ok = []
    for sp in sorted(species):
        Ls = {s: v for s, v in L.items() if s.startswith(sp + '|')}
        Rs = {s: v for s, v in R.items() if s.startswith(sp + '|')}
        if not Ls or not Rs:
            continue
        lb = max(Ls, key=lambda s: Ls[s][2])
        rb = max(Rs, key=lambda s: Rs[s][2])
        cov_l = (abs(Ls[lb][1] - Ls[lb][0]) + 1) / Ls[lb][3]
        cov_r = (abs(Rs[rb][1] - Rs[rb][0]) + 1) / Rs[rb][3]
        if lb != rb and lb not in Rs and rb not in Ls and cov_l >= MIN_COV and cov_r >= MIN_COV:
            ok.append(f'{sp}:{lb.split("|", 1)[1]}({cov_l:.2f})/{rb.split("|", 1)[1]}({cov_r:.2f})')
    return ok


rows = [l.rstrip('\n').split('\t') for l in open(junc_path)]
head = rows[0]
ix = {k: i for i, k in enumerate(head)}
out = open(out_path, 'w')
out.write('\t'.join(head + ['class', 'tag', 'split_species', 'n_split_species', 'split']) + '\n')
cnt = collections.defaultdict(collections.Counter)
tagged = collections.defaultdict(set)
n_split = 0
split_tx = set()
for r in rows[1:]:
    la, ra, sup = int(r[ix['left_aa']]), int(r[ix['right_aa']]), r[ix['support']]
    c = 'segment_lt_30aa' if la < 30 or ra < 30 else classify(r[ix['left_seg']], r[ix['right_seg']])
    cnt[sup][c] += 1
    tag = ''
    if c == 'different_proteins' and sup == 'unsupported':
        tag = 'possible_gap_fusion'
    elif c == 'same_protein_noncollinear' and sup == 'unsupported':
        tag = 'possible_gap_misorder'
    sp = split_species(r[ix['left_seg']], r[ix['right_seg']]) if (sup == 'unsupported' and c != 'segment_lt_30aa') else []
    split = 'yes' if len(sp) >= MIN_SPECIES else 'no'
    if split == 'yes':
        n_split += 1
        split_tx.add(r[ix['transcript']])
    if tag:
        tagged[tag].add(r[ix['transcript']].rsplit('.', 1)[0])
    out.write('\t'.join(r + [c, tag, ';'.join(sp), str(len(sp)), split]) + '\n')
out.close()

order = ['same_protein_collinear', 'same_protein_noncollinear', 'different_proteins',
         'one_side_only', 'no_hit', 'segment_lt_30aa']
print(f"{'class':28s}" + ''.join(f'{k:>14s}' for k in ('rna', 'protein_only', 'unsupported')))
for c in order:
    print(f'{c:28s}' + ''.join(f'{cnt[k][c]:14d}' for k in ('rna', 'protein_only', 'unsupported')))
for t in ('possible_gap_fusion', 'possible_gap_misorder'):
    n_j = sum(1 for r in open(out_path) if f'\t{t}\t' in r)
    print(f'{t}: {n_j} junctions in {len(tagged[t])} genes')
print(f'split (strict rule, >= {MIN_SPECIES} proteomes, cov >= {MIN_COV}): {n_split} junctions in '
      f'{len(split_tx)} transcripts, {len({t.rsplit(".", 1)[0] for t in split_tx})} genes')
