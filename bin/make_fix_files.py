#!/usr/bin/env python3
"""Turn the reviewed product curation (refined gene IDs) into funannotate's
--fix and --remove files keyed by the post-rename locus tags.

The locus tags are computed with funannotate's own convertgff2tbl on the same GFF3 and FASTA that
pass 2 uses, with the same prefix, so they are exactly the tags pass 2 assigns. Run inside the
funannotate 1.8.17 image.
usage: make_fix_files.py refined.gff3 scaffolds.fa LOCUS_PREFIX product_curation.tsv fix.tsv remove.txt map.tsv
"""
import os
import sys
import tempfile

import logging

from funannotate import library as lib

logging.basicConfig(level=logging.WARNING)
lib.log = logging.getLogger('funannotate')   # annotate normally sets this up itself

gff, fasta, prefix, cur_p, fix_p, rem_p, map_p = sys.argv[1:8]
# as funannotate annotate does with --rename
prefix = prefix.replace('_', '')
tmp = tempfile.mkdtemp(dir='.')
_, geneDB = lib.convertgff2tbl(gff, prefix, fasta, os.path.join(tmp, 'p.fa'), os.path.join(tmp, 't.fa'),
                               os.path.join(tmp, 'g.tbl'), external=True)
gene2locus = {}
for tx, locus in geneDB.items():
    gene2locus.setdefault(tx.rsplit('.', 1)[0], locus)
with open(map_p, 'w') as o:
    o.write('refined_gene\tlocus_tag\n')
    for g, l in sorted(gene2locus.items(), key=lambda x: x[1]):
        o.write(f'{g}\t{l}\n')
n_fix = n_rem = missing = 0
with open(cur_p) as fi, open(fix_p, 'w') as ff, open(rem_p, 'w') as fr:
    for line in fi:
        f = line.rstrip('\n').split('\t')
        if f[0] == 'gene' or line.startswith('#'):
            continue
        g, sym, new, action = f[0], f[1], f[3], f[4]
        locus = gene2locus.get(g)
        if not locus:
            missing += 1
            continue
        if action == 'fix':
            ff.write(f'{locus}\t{sym}\t{new}\n')
            n_fix += 1
        else:
            fr.write(f'{locus}\n')
            n_rem += 1
print(f'fix {n_fix}, remove {n_rem}, genes not in the final GFF3 (removed by screening) {missing}')
