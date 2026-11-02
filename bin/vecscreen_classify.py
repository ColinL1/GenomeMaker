#!/usr/bin/env python3
"""Classify VecScreen (blastn vs UniVec_Core) hits with NCBI's VecScreen match categories.

Raw blastn score thresholds (VecScreen with its standard parameters):
  terminal match (within 25 bp of a sequence end): Strong >= 24, Moderate 19-23, Weak 16-18
  internal match:                                  Strong >= 30, Moderate 25-29, Weak 23-24
Overlapping hits on a query are merged and the region takes its best category.

usage: vecscreen_classify.py hits.tsv out.tsv
hits.tsv: blastn -outfmt "6 qseqid qstart qend sseqid sstart send pident length evalue bitscore score qlen"
Prints counts per category; out.tsv lists every classified region.
"""
import collections
import sys

TERM = 25
RANK = {'Strong': 3, 'Moderate': 2, 'Weak': 1}


def category(score, terminal):
    if terminal:
        return 'Strong' if score >= 24 else 'Moderate' if score >= 19 else 'Weak' if score >= 16 else None
    return 'Strong' if score >= 30 else 'Moderate' if score >= 25 else 'Weak' if score >= 23 else None


regions = collections.defaultdict(list)
for line in open(sys.argv[1]):
    f = line.rstrip('\n').split('\t')
    q, a, b, s, score, qlen = f[0], int(f[1]), int(f[2]), f[3], int(f[10]), int(f[11])
    a, b = min(a, b), max(a, b)
    terminal = a <= TERM or b > qlen - TERM
    c = category(score, terminal)
    if c:
        regions[q].append([a, b, c, terminal, s, score])
cnt = collections.Counter()
with open(sys.argv[2], 'w') as o:
    o.write('seq\tstart\tend\tcategory\tterminal\tbest_subject\tbest_score\n')
    for q, lst in regions.items():
        lst.sort()
        merged = []
        for r in lst:
            if merged and r[0] <= merged[-1][1]:
                m = merged[-1]
                m[1] = max(m[1], r[1])
                if RANK[r[2]] > RANK[m[2]] or (RANK[r[2]] == RANK[m[2]] and r[5] > m[5]):
                    m[2], m[3], m[4], m[5] = r[2], r[3], r[4], r[5]
            else:
                merged.append(list(r))
        for a, b, c, t, s, sc in merged:
            cnt[(c, 'terminal' if t else 'internal')] += 1
            o.write(f'{q}\t{a}\t{b}\t{c}\t{int(t)}\t{s}\t{sc}\n')
for k in sorted(cnt):
    print(f'{k[0]}\t{k[1]}\t{cnt[k]}')
if not cnt:
    print('no VecScreen matches')
