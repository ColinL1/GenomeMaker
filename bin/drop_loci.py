#!/usr/bin/env python3
"""Remove genes from funannotate's GFF3, annotations.txt or GenBank file, so the release files match the
NCBI submission (fix_tbl.py removes rRNA fragments from the .tbl and lists them as action 'drop' in rrna.tsv).

The format comes from the extension: .gff3 (gene ID and every feature whose ID/Parent is the gene or one
of its transcripts <gene>-T<n>), .annotations.txt (GeneID column), .gbk (feature blocks with that
/locus_tag). Lines outside these features are copied unchanged.

usage: drop_loci.py rrna.tsv in_file out_file
"""
import re
import sys

tsv, fin, fout = sys.argv[1:4]
drop = set()
with open(tsv) as fh:
    head = fh.readline().rstrip('\n').split('\t')
    if 'locus_tag' in head:
        for line in fh:
            r = dict(zip(head, line.rstrip('\n').split('\t')))
            if r.get('action') == 'drop':
                drop.add(r['locus_tag'])

n = 0
with open(fin) as fi, open(fout, 'w') as o:
    if fin.endswith('.gff3'):
        for line in fi:
            f = line.split('\t')
            if len(f) > 8:
                ids = re.findall(r'(?:^|;)(?:ID|Parent)=([^;\n]+)', f[8])
                if any(re.sub(r'-T\d+$', '', x) in drop for i in ids for x in i.split(',')):
                    n += f[2] == 'gene'
                    continue
            o.write(line)
    elif fin.endswith('.annotations.txt'):
        for line in fi:
            if line.split('\t', 1)[0] in drop:
                n += 1
                continue
            o.write(line)
    elif fin.endswith('.gbk'):
        block = []          # current feature block (key line + qualifier lines)

        def flush():
            global n
            m = re.search(r'/locus_tag="([^"]+)"', ''.join(block))
            if m and m.group(1) in drop:
                n += block[0][5:21].strip() == 'gene'
            else:
                o.writelines(block)
            block.clear()
        for line in fi:
            if line.startswith('     ') and line[5] != ' ':      # feature key line
                flush()
                block.append(line)
            elif line.startswith(' ' * 21) and block:          # qualifier or location continuation
                block.append(line)
            else:
                flush()
                o.write(line)
        flush()
    else:
        sys.exit(f'ERROR: unknown format {fin}')
print(f'{fin}: {n} of {len(drop)} genes removed')
if n != len(drop):
    sys.exit(f'ERROR: {fin}: expected to remove {len(drop)} genes, removed {n}')
