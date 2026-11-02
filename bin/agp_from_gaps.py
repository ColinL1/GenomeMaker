#!/usr/bin/env python3
"""AGP of gap-filled scaffolds, from their N runs.

TGS-GapCloser and Medaka close some of RagTag's gaps, so RagTag's AGP no longer matches the
sequence. Every N run left becomes a gap line; the sequence between N runs a W line.
Gap lines keep RagTag's type and linkage columns (U / 'scaffold yes align_genus' by default,
N with ragtag_infer_gaps), read from RagTag's AGP; the gap length is the N-run length.
Component ids are <scaffold>_c<k> (the gap-filled pieces are no longer RagTag's contigs).

usage: agp_from_gaps.py scaffolds.fasta ragtag.scaffold.agp out.agp
"""
import re
import sys

fa, ragtag_agp, out_p = sys.argv[1:4]


def records(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith('>'):
            if name is not None:
                yield name, ''.join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        yield name, ''.join(buf)


# RagTag's gap columns: 5 type (U/N), 7 gap type, 8 linkage, 9 linkage evidence
styles = set()
for l in open(ragtag_agp):
    x = l.rstrip('\n').split('\t')
    if not l.startswith('#') and len(x) >= 9 and x[4] in ('N', 'U'):
        styles.add((x[4], x[6], x[7], x[8]))
if len(styles) > 1:
    sys.exit(f'ERROR: RagTag AGP has gap lines of more than one kind: {sorted(styles)}')
gap_type, gap_kind, linkage, evidence = styles.pop() if styles else (None, None, None, None)

n_gaps = 0
with open(out_p, 'w') as o:
    o.write('## agp-version 2.1\n# AGP rebuilt from the N runs of the gap-filled scaffolds (agp_from_gaps.py)\n')
    for name, seq in records(fa):
        part = comp = 0
        pos = 0
        for m in list(re.finditer(r'[Nn]+', seq)) + [None]:
            start, end = (m.start(), m.end()) if m else (len(seq), len(seq))
            if start > pos:
                part += 1
                comp += 1
                o.write(f'{name}\t{pos + 1}\t{start}\t{part}\tW\t{name}_c{comp}\t1\t{start - pos}\t+\n')
            if m:
                if gap_type is None:
                    sys.exit(f'ERROR: {name} has an N run at {start + 1}, but the RagTag AGP has no gap lines to take the gap type from')
                part += 1
                n_gaps += 1
                o.write(f'{name}\t{start + 1}\t{end}\t{part}\t{gap_type}\t{end - start}\t{gap_kind}\t{linkage}\t{evidence}\n')
            pos = end
print(f'{n_gaps} gaps')
