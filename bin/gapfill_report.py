#!/usr/bin/env python3
"""Gap filling (03e) in numbers, as a MultiQC custom-content table (gapfill_mqc.tsv).

Rows: RagTag scaffolds -> TGS-GapCloser -> Medaka. Gaps and contig N50 from gfastats; joins and
fills from TGS-GapCloser's <prefix>.gap_fill_detail, where a scaffold is a list of segments
  start end S orig_start orig_end   sequence of the input scaffold (orig_*: input coordinates)
  start end F                       fill (read sequence)
  start end N                       gap left open
Two S segments in a row are a direct join (the contig ends overlapped); the bases dropped there
are the input coordinates between the two segments minus RagTag's gap (from its AGP).
Medaka closes the N runs that reads span: its gap count is the gfastats one.

usage: gapfill_report.py gap_fill_detail ragtag.scaffold.agp ragtag.gfastats.txt tgs.gfastats.txt
                         medaka.gfastats.txt|NO_* gapfill_mqc.tsv
"""
import collections
import os
import sys

detail, agp, gf_ragtag, gf_tgs, gf_medaka, out_p = sys.argv[1:7]


def gfastats(path):
    if os.path.basename(path).startswith('NO_'):
        return None
    d = {}
    for l in open(path):
        if ':' in l:
            k, v = l.split(':', 1)
            d[k.strip()] = v.strip()
    return d


ragtag_gaps = collections.defaultdict(list)   # object -> [(start, end)] of RagTag gaps
for l in open(agp):
    x = l.rstrip('\n').split('\t')
    if not l.startswith('#') and len(x) >= 6 and x[4] in ('N', 'U'):
        ragtag_gaps[x[0]].append((int(x[1]), int(x[2])))

segs = collections.defaultdict(list)
name = None
for l in open(detail):
    if l.startswith('>'):
        name = l[1:].strip()
        continue
    f = l.split()
    if len(f) >= 3:
        segs[name].append(f)

joins = trimmed = fills = fill_bp = 0
for name, ss in segs.items():
    for prev, cur in zip(ss, ss[1:]):
        if cur[2] == 'F':
            fills += 1
            fill_bp += int(cur[1]) - int(cur[0]) + 1
        elif cur[2] == 'S' and prev[2] == 'S':
            joins += 1
            lo, hi = int(prev[4]), int(cur[3])
            gap_bp = sum(min(e, hi - 1) - max(s, lo + 1) + 1 for s, e in ragtag_gaps[name] if s < hi and e > lo)
            trimmed += hi - lo - 1 - gap_bp

stats = [('RagTag scaffolds', gfastats(gf_ragtag)), ('TGS-GapCloser', gfastats(gf_tgs)),
         ('Medaka', gfastats(gf_medaka))]
with open(out_p, 'w') as o:
    o.write("# id: 'gapfill'\n")
    o.write("# section_name: '03e Gap filling'\n")
    o.write('# description: "Gaps of the RagTag scaffolds closed by TGS-GapCloser (--ne, ONT reads) and by the '
            'second Medaka pass, and the contig N50 after each step. Direct joins: contig ends overlapped and were '
            'joined without fill sequence; overlap bp trimmed: bases removed there. Assembly that goes on: '
            f'{"Medaka" if stats[2][1] else "TGS-GapCloser"} output."\n')
    o.write("# plot_type: 'table'\n")
    o.write("# pconfig:\n#     id: 'gapfill_table'\n#     col1_header: 'Stage'\n#     sort_rows: false\n")
    o.write('Stage\tGaps\tContig N50\tDirect joins\tOverlap bp trimmed\tFills\tFill bp\n')
    for stage, g in stats:
        tgs = [joins, trimmed, fills, fill_bp] if stage == 'TGS-GapCloser' else ['', '', '', '']
        if g is None:
            o.write(f'{stage}\tnot run\tnot run\t\t\t\t\n')
            continue
        o.write('\t'.join(str(v) for v in [stage, g.get('# gaps in scaffolds', 'NA'), g.get('Contig N50', 'NA')] + tgs) + '\n')
print(f'direct joins {joins} ({trimmed} bp trimmed); fills {fills} ({fill_bp} bp)')
