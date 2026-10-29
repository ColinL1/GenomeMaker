#!/usr/bin/env python3
"""Release GFF3 and annotations.txt made to agree with the final .tbl (the one table2asn submits).

FIX_TBL corrects funannotate's .tbl (gene ends run out to gaps, product clean-ups and curation, EC
numbers, gene symbols); funannotate's GFF3 and annotations.txt are written before that. This copies,
per gene, from the final .tbl:
  - GFF3 gene: start/end (partial marks dropped), Name (the gene symbol; removed when the .tbl has none,
    protein-coding genes only);
  - GFF3 mRNA: product (the mRNA product), EC_number (the CDS's EC numbers, ',' joined);
  - annotations.txt mRNA rows: Name, Product, EC_number (';' joined).
GFF3 mRNAs are matched by ID (the .tbl transcript_id without 'gnl|...|' and '_mrna'); annotations.txt
rows by gene and funannotate's ', variant N' product suffix: annotations.txt numbers isoforms in its own
order and its rows carry the gene span, so the suffix is the only link. Isoforms without a suffix
(e.g. several 'hypothetical protein') share a key and must carry identical values in the .tbl; the
number of rows per key must equal the number of .tbl mRNAs. Genes in one file but not the other stop the run (rRNA fragments are removed from the GFF3 and annotations.txt by drop_loci.py before).

usage: sync_release_files.py final.tbl in.gff3 in.annotations.txt out.gff3 out.annotations.txt report.tsv
"""
import collections
import re
import sys

tbl_p, gff_in, ann_in, gff_out, ann_out, rep = sys.argv[1:7]


def esc(v):
    return v.replace('%', '%25').replace(';', '%3B').replace('=', '%3D').replace('&', '%26')


# ---- final .tbl: genes (locus_tag -> span, symbol, mRNA list) and mRNAs (ID -> product, ECs, span)
genes = collections.OrderedDict()
mrna = {}
cur_gene = cur = None
for line in open(tbl_p):
    line = line.rstrip('\n')
    if line.startswith('>Feature'):
        cur_gene = cur = None
        continue
    f = line.split('\t')
    if line and not line.startswith('\t'):
        a, b = int(f[0].lstrip('<>')), int(f[1].lstrip('<>'))
        if len(f) >= 3 and f[2]:
            cur = {'type': f[2], 'span': [min(a, b), max(a, b)], 'q': collections.defaultdict(list)}
            if f[2] == 'gene':
                cur_gene = cur
                cur['kids'] = []
            elif cur_gene is not None:
                cur_gene['kids'].append(cur)
        elif cur is not None:                       # further intervals of a joined location
            cur['span'] = [min(cur['span'][0], a, b), max(cur['span'][1], a, b)]
        continue
    if cur is not None and line.startswith('\t\t\t') and len(f) > 4:
        cur['q'][f[3]].append(f[4])
        if cur is cur_gene and f[3] == 'locus_tag':
            genes[f[4]] = cur

for lt, g in genes.items():
    g['mrnas'] = []
    for k in g['kids']:
        tid = re.sub(r'^gnl\|[^|]+\|', '', (k['q'].get('transcript_id') or [''])[0]).rsplit('_mrna', 1)[0]
        if not tid:
            continue
        if k['type'] == 'mRNA':
            mrna[tid] = {'product': (k['q'].get('product') or [''])[0], 'ec': []}
            g['mrnas'].append(tid)
        elif k['type'] == 'CDS' and tid in mrna:
            mrna[tid]['ec'] = k['q'].get('EC_number', [])

cnt = collections.Counter()


def set_attr(attrs, key, val):
    """attrs: list of 'k=v'; val None removes the key. Returns True when changed."""
    old = next((x for x in attrs if x.split('=', 1)[0] == key), None)
    new = None if val in (None, '') else f'{key}={esc(val)}'
    if old == new or (old and new and old.split('=', 1)[1] == val):
        return False
    if old is not None:
        i = attrs.index(old)
        if new is None:
            del attrs[i]
        else:
            attrs[i] = new
    elif new is not None:
        attrs.append(new)
    return True


# ---- GFF3
seen = set()
out = []
for line in open(gff_in):
    f = line.rstrip('\n').split('\t')
    if line.startswith('#') or len(f) < 9:
        out.append(line.rstrip('\n'))
        continue
    attrs = [x for x in f[8].rstrip(';').split(';') if x]
    ida = next((x.split('=', 1)[1] for x in attrs if x.startswith('ID=')), '')
    if f[2] == 'gene':
        g = genes.get(ida)
        if g is None:
            sys.exit(f'ERROR: gene {ida} of the GFF3 is not in the final .tbl')
        seen.add(ida)
        if [int(f[3]), int(f[4])] != g['span']:
            f[3], f[4] = map(str, g['span'])
            cnt['gff_gene_span_changed'] += 1
        if g['mrnas']:
            cnt['gff_gene_name_changed'] += set_attr(attrs, 'Name', (g['q'].get('gene') or [None])[0])
    elif f[2] == 'mRNA' and ida in mrna:
        m = mrna[ida]
        cnt['gff_mrna_product_changed'] += set_attr(attrs, 'product', m['product'])
        cnt['gff_mrna_EC_changed'] += set_attr(attrs, 'EC_number', ','.join(m['ec']))
    elif f[2] == 'mRNA':
        sys.exit(f'ERROR: mRNA {ida} of the GFF3 is not in the final .tbl')
    f[8] = ';'.join(attrs) + ';'
    out.append('\t'.join(f))
missing = [lt for lt in genes if lt not in seen]
if missing:
    sys.exit(f'ERROR: {len(missing)} genes of the final .tbl are not in the GFF3: {", ".join(missing[:5])}')
with open(gff_out, 'w') as o:
    o.write('\n'.join(out) + '\n')

# ---- annotations.txt (mRNA rows; isoforms in the same order as the .tbl's mRNAs)
rows = [l.rstrip('\n').split('\t') for l in open(ann_in)]
hd = {k: i for i, k in enumerate(rows[0])}


def variant(product):
    v = re.search(r', variant (\d+)$', product)
    return v.group(1) if v else ''


by_variant = collections.defaultdict(list)
for lt, g in genes.items():
    for tid in g['mrnas']:
        by_variant[(lt, variant(mrna[tid]['product']))].append(mrna[tid])
for key, ms in by_variant.items():
    if len({(m['product'], tuple(m['ec'])) for m in ms}) > 1:
        sys.exit(f'ERROR: gene {key[0]}: isoforms without a distinguishing variant number differ in product or EC')
used = collections.Counter()
for r in rows[1:]:
    if r[hd['Feature']] != 'mRNA':
        continue
    key = (r[hd['GeneID']], variant(r[hd['Product']]))
    used[key] += 1
    if key not in by_variant or used[key] > len(by_variant[key]):
        sys.exit(f"ERROR: annotations.txt row {r[hd['TranscriptID']]} has no matching mRNA in the final .tbl")
    g, m = genes[key[0]], by_variant[key][0]
    for col, val in (('Name', (g['q'].get('gene') or [''])[0]), ('Product', m['product']),
                     ('EC_number', ';'.join(m['ec']))):
        if r[hd[col]] != val:
            r[hd[col]] = val
            cnt[f'annotations_{col}_changed'] += 1
short = [k for k, ms in by_variant.items() if used[k] != len(ms)]
if short:
    sys.exit(f'ERROR: {len(short)} genes have mRNAs of the final .tbl without a row in annotations.txt: '
             f'{", ".join(k[0] for k in short[:5])}')
with open(ann_out, 'w') as o:
    o.write('\n'.join('\t'.join(r) for r in rows) + '\n')

with open(rep, 'w') as o:
    for k in sorted(cnt):
        o.write(f'{k}\t{cnt[k]}\n')
print(dict(cnt))
