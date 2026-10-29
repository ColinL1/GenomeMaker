#!/usr/bin/env python3
"""Add Class=<class/family> to rmOutToGFF3.pl output, from the (relabelled) RepeatMasker .out.

Features and .out records are matched in order; the family in Target= must equal the .out family,
and the counts must be equal, otherwise the script fails.
usage: repeat_gff_classes.py repeatmasker.out in.gff3 out.gff3
"""
import re
import sys

recs = []
for line in open(sys.argv[1]):
    f = line.split()
    if len(f) >= 11 and f[0].isdigit():
        recs.append((f[9], f[10]))
n = 0
with open(sys.argv[3], 'w') as o:
    for line in open(sys.argv[2]):
        if line.startswith('#') or not line.strip():
            o.write(line)
            continue
        fam = re.search(r'Target=(\S+)', line).group(1)
        if n >= len(recs) or recs[n][0] != fam:
            sys.exit(f'ERROR: feature {n + 1} family {fam} does not match .out record {recs[n] if n < len(recs) else None}')
        o.write(line.rstrip('\n') + f';Class={recs[n][1]}\n')
        n += 1
if n != len(recs):
    sys.exit(f'ERROR: {n} GFF features but {len(recs)} .out records')
print(f'{n} features = {len(recs)} .out records')
