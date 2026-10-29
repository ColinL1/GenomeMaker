#!/usr/bin/env python3
"""Mechanical NCBI fixes to funannotate's .tbl before table2asn.

  - ncRNA features get ncRNA_class from the refined GFF3 (matched by scaffold, start, end, strand);
  - EC_number qualifiers are removed from CDS/proteins whose product is 'hypothetical protein'
    (SEQ_FEAT.BadProteinName: unknown or hypothetical protein should not have an EC number);
  - incomplete EC numbers are padded with '-' to four fields, '2.7.7' -> '2.7.7.-'
    (SEQ_FEAT.BadEcNumberFormat);
  - trailing '.', ',', ';' or ':' are removed from product names (SEQ_FEAT.BadTrailingCharacter);
  - rRNA products 'Eukaryotic large/small subunit ribosomal RNA' become '28S/18S ribosomal RNA';
    rRNA models shorter than 90% of the Rfam model (LSU RF02543 3401, SSU RF01960 1831, 5.8S RF00002
    154, 5S RF00001 119 bp) are fragments. NCBI accepts a partial
    end only at a sequence end or a gap, so per fragment, on the rRNA and on its gene:
      an end at a sequence end or next to a gap is marked partial ('<' 5', '>' 3'), the other end is not;
      a fragment with no such end is removed with its gene (unmarked it fails the discrepancy report,
      FATAL SHORT_RRNA; marked it fails validation, PartialProblem);
    every rRNA and the decision is listed in rrna.tsv (release_sets.py removes the same genes);
  - protein-coding genes span exactly their mRNAs; a gene end is partial only if every mRNA reaching it is
    partial there (NCBI GENE_PARTIAL_CONFLICT);
  - CDS products lose funannotate's ', variant N' (isoforms share the protein name; it stays on the mRNA);
  - EC numbers table2asn changed in a preflight run (its .ecn: 'replaced' by one number -> that number;
    'split' or 'deleted' -> removed), so the final table2asn run changes none;
  - products from the final product curation (product_curation.tsv, keyed by gene) are applied to the
    gene's mRNAs and CDSs ('remove' also drops the gene symbol); a product table2asn still rewrote in the
    preflight run (its .fixedproducts) gets table2asn's own name, so the final run rewrites none;
  - human chromosome-ORF gene symbols (CXorf38 ...) and symbols containing 'like' (PC3LIKE) are removed;
  - mRNA and CDS features of genes whose span equals another gene's span get gene cross-references
    (locus_tag, and gene symbol if any), so they are unambiguous (SEQ_FEAT.MissingGeneXref);
  - with --strip-gene-suffix: funannotate's '_N' suffix (added when one symbol is assigned to several
    genes) is removed from every gene symbol, so /gene is the symbol itself (NCBI BAD_GENE_NAME: > 10
    characters only because of the suffix); locus_tags stay the unique IDs.
Nothing else is changed. Every change is counted in the report.

usage: fix_tbl.py [--strip-gene-suffix] in.tbl refined.gff3 out.tbl report.tsv scaffolds.fa rrna.tsv
       [x.ecn product_curation.tsv x.fixedproducts]
(the last three come from the preflight table2asn run; NO_* placeholders or missing files mean none)
"""
import collections
import os
import re
import sys

STRIP_SUFFIX = '--strip-gene-suffix' in sys.argv
argv = [x for x in sys.argv if x != '--strip-gene-suffix']
tbl_in, gff, tbl_out, rep, fasta, rrna_out = argv[1:7]


def opt(i):
    p = argv[i] if len(argv) > i else ''
    return p if p and not os.path.basename(p).startswith('NO_') and os.path.exists(p) else None


ecn_p, cur_p, fxp_p = opt(7), opt(8), opt(9)

ncls = {}
tx_gene = {}                     # refined mRNA -> gene
for line in open(gff):
    if line.startswith('#'):
        continue
    f = line.rstrip('\n').split('\t')
    if len(f) > 8 and f[2] == 'ncRNA':
        m = re.search(r'ncRNA_class=([^;]+)', f[8])
        if m:
            ncls[(f[0], int(f[3]), int(f[4]), f[6])] = m.group(1)
    elif len(f) > 8 and f[2] == 'mRNA':
        tx_gene[re.search(r'ID=([^;]+)', f[8]).group(1)] = re.search(r'Parent=([^;]+)', f[8]).group(1)

cnt = collections.Counter()
out = []
lines = open(tbl_in).read().split('\n')
seq = None
i = 0
# features: a location line (no leading tab) followed by qualifier lines ('\t\t\t')
blocks = []
cur = None
for line in lines:
    if line.startswith('>Feature'):
        seq = line.split()[1]
        if cur:
            blocks.append(cur)
            cur = None
        blocks.append(('header', [line]))
        continue
    if line and not line.startswith('\t'):
        f = line.split('\t')
        if len(f) >= 3 and f[2]:          # new feature
            if cur:
                blocks.append(cur)
            cur = ('feat', [line], seq)
            continue
    if cur:
        cur[1].append(line)
    else:
        blocks.append(('other', [line]))
if cur:
    blocks.append(cur)

if STRIP_SUFFIX:
    for b in blocks:
        for k, l in enumerate(b[1]):
            if l.startswith('\t\t\tgene\t') and re.search(r'_\d+$', l):
                b[1][k] = re.sub(r'_\d+$', '', l)
                cnt['gene_suffix_removed'] += 1

RRNA = {'Eukaryotic large subunit ribosomal RNA': ('28S ribosomal RNA', 3401),
        'Eukaryotic small subunit ribosomal RNA': ('18S ribosomal RNA', 1831),
        '5.8S ribosomal RNA': ('5.8S ribosomal RNA', 154),     # Rfam RF00002
        '5S ribosomal RNA': ('5S ribosomal RNA', 119)}         # Rfam RF00001


def set_partial(loc_line, five, three):
    f = loc_line.split('\t')
    f[0], f[1] = ('<' if five else '') + f[0].lstrip('<>'), ('>' if three else '') + f[1].lstrip('<>')
    return '\t'.join(f)


def loc(b):
    f = b[1][0].split('\t')
    return int(f[0].lstrip('<>')), int(f[1].lstrip('<>'))


# sequences of the scaffolds that carry rRNA fragments (sequence ends and gaps)
need = set()
for b in blocks:
    if b[0] == 'feat' and b[1][0].split('\t')[2] == 'rRNA':
        need.add(b[2])
seqs, name = {}, None
for line in open(fasta):
    if line.startswith('>'):
        name = line[1:].split()[0]
        if name in need:
            seqs[name] = []
    elif name in seqs:
        seqs[name].append(line.strip())
seqs = {k: ''.join(v) for k, v in seqs.items()}


def at_end_or_gap(sid, pos):
    """pos is the base just outside the feature end."""
    s = seqs[sid]
    return pos < 1 or pos > len(s) or s[pos - 1] in 'Nn'


# rRNA renaming and partial marks are applied to the block list before the main pass
rrna_rows, drop = [], set()
last_gene = None
for i, b in enumerate(blocks):
    if b[0] != 'feat':
        continue
    ftype = b[1][0].split('\t')[2]
    if ftype == 'gene':
        last_gene = i
    if ftype != 'rRNA':
        continue
    prod = next((l.split('\t')[4] for l in b[1] if l.startswith('\t\t\tproduct\t')), None)
    if prod not in RRNA:
        continue
    name, model_len = RRNA[prod]
    b[1][:] = [f'\t\t\tproduct\t{name}' if l.startswith('\t\t\tproduct\t') else l for l in b[1]]
    cnt['rRNA_renamed'] += name != prod
    p5, p3 = loc(b)
    length = abs(p3 - p5) + 1
    gb = blocks[last_gene][1] if last_gene is not None else None
    if gb is None or sorted(loc(blocks[last_gene])) != sorted((p5, p3)):
        sys.exit(f'ERROR: rRNA {b[2]}:{p5}-{p3} has no gene with the same span before it')
    tag = next((l.split('\t')[4] for l in gb if l.startswith('\t\t\tlocus_tag\t')), '')
    step = 1 if p5 <= p3 else -1
    five = three = False
    if length >= 0.9 * model_len or len([l for l in b[1] if l and not l.startswith('\t')]) != 1:
        action = 'complete'
    else:
        five, three = at_end_or_gap(b[2], p5 - step), at_end_or_gap(b[2], p3 + step)
        if five or three:
            action = 'partial_at_sequence_end_or_gap'
        else:
            action = 'drop'
            drop |= {i, last_gene}
    if action != 'drop':
        b[1][0] = set_partial(b[1][0], five, three)
        gb[0] = set_partial(gb[0], five, three)
    cnt[f'rRNA_{action}'] += 1
    rrna_rows.append((tag, b[2], min(p5, p3), max(p5, p3), '+' if step == 1 else '-', name, length, model_len,
                      f'{100 * length / model_len:.0f}', '<' if five else '', '>' if three else '', action))
blocks = [b for i, b in enumerate(blocks) if i not in drop]
with open(rrna_out, 'w') as o:
    o.write('locus_tag\tseqid\tstart\tend\tstrand\tproduct\tlength\tmodel_length\tpct_model\tpartial_5\tpartial_3\taction\n')
    for r in rrna_rows:
        o.write('\t'.join(map(str, r)) + '\n')

# protein-coding genes span exactly their mRNAs, and a gene end is partial only if every mRNA reaching that
# end is partial there (a partial gene end beside a complete mRNA is NCBI GENE_PARTIAL_CONFLICT)
groups = []                      # [gene block index, [child block indexes]]
for i, b in enumerate(blocks):
    if b[0] != 'feat':
        continue
    if b[1][0].split('\t')[2] == 'gene':
        groups.append([i, []])
    elif groups:
        groups[-1][1].append(i)
for gi, kids in groups:
    ends = []                    # per mRNA: (5' position, 5' partial, 3' position, 3' partial)
    for k in kids:
        ls = blocks[k][1]
        if ls[0].split('\t')[2] != 'mRNA':
            continue
        locs = [l.split('\t') for l in ls if l and not l.startswith('\t')]
        ends.append((int(locs[0][0].lstrip('<>')), locs[0][0].startswith('<'),
                     int(locs[-1][1].lstrip('<>')), locs[-1][1].startswith('>')))
    if not ends:
        continue
    f = blocks[gi][1][0].split('\t')
    g5, g3 = int(f[0].lstrip('<>')), int(f[1].lstrip('<>'))
    plus = ends[0][0] <= ends[0][2]
    n5 = min(e[0] for e in ends) if plus else max(e[0] for e in ends)
    n3 = max(e[2] for e in ends) if plus else min(e[2] for e in ends)
    p5 = all(e[1] for e in ends if e[0] == n5)
    p3 = all(e[3] for e in ends if e[2] == n3)
    new = [('<' if p5 else '') + str(n5), ('>' if p3 else '') + str(n3)]
    if new != f[:2]:
        if (n5, n3) != (g5, g3):
            cnt['gene_span_set_to_mrnas'] += 1
        if (p5, p3) != (f[0].startswith('<'), f[1].startswith('>')):
            cnt['gene_partial_set_from_mrnas'] += 1
        f[0], f[1] = new
        blocks[gi][1][0] = '\t'.join(f)

# products from the final curation (by gene) and table2asn's own rewrites from the preflight run (by protein)
curation = {}
for l in (open(cur_p) if cur_p else []):
    f = l.rstrip('\n').split('\t')
    if len(f) >= 5 and f[0] != 'gene':
        curation[f[0]] = (f[3], f[4])
rewrites = {}
for l in (open(fxp_p) if fxp_p else []):
    m = re.match(r'Changed (.+) gnl\|[^|]+\|([^:\s]+)', l.strip())
    if m:
        rewrites[m.group(2)] = m.group(1)       # 'OLD to NEW'; split once the current product is known


def qual(ls, key):
    return next((l.split('\t')[4] for l in ls if l.startswith(f'\t\t\t{key}\t')), None)


def set_qual(ls, key, val):
    ls[:] = [f'\t\t\t{key}\t{val}' if l.startswith(f'\t\t\t{key}\t') else l for l in ls]


for gi, kids in groups:
    for k in kids:
        ls = blocks[k][1]
        ftype = ls[0].split('\t')[2]
        if ftype not in ('mRNA', 'CDS'):
            continue
        tid = re.sub(r'^gnl\|[^|]+\|', '', qual(ls, 'transcript_id') or qual(ls, 'protein_id') or '').rsplit('_mrna', 1)[0]
        prod = qual(ls, 'product')
        if prod is None or not tid:
            continue
        variant = re.search(r', variant \d+$', prod) if ftype == 'mRNA' else None
        base = prod[:variant.start()] if variant else prod
        new = base
        hit = curation.get(tx_gene.get(tid, tid.rsplit('.t', 1)[0]))
        if hit:
            new = hit[0]
            if hit[1] == 'remove' and qual(blocks[gi][1], 'gene') is not None:
                blocks[gi][1][:] = [l for l in blocks[gi][1] if not l.startswith('\t\t\tgene\t')]
                cnt['curation_symbol_removed'] += 1
        rw = rewrites.get(tid)
        if rw and rw.startswith(f'{new} to '):
            new = rw[len(new) + 4:]
            cnt[f'table2asn_rewrite_applied_{ftype}'] += 1
        elif rw and ftype == 'CDS':
            cnt['table2asn_rewrite_not_matched'] += 1
        if new != base:
            set_qual(ls, 'product', new + (variant.group(0) if variant else ''))
            cnt[f'curation_product_changed_{ftype}'] += 1

# EC numbers table2asn changed in the preflight run (.ecn: protein, 'EC number replaced/split/deleted',
# old number, new numbers): one new number replaces the old one, a split or a deletion removes it
EC_CHANGE = {}
for l in (open(ecn_p) if ecn_p else []):
    f = l.rstrip('\n').split('\t')
    if len(f) >= 3 and f[2]:
        new = [x for x in f[3:] if x.strip()]
        EC_CHANGE[f[2]] = new[0] if len(new) == 1 else None
# Human chromosome-ORF symbols (C1orf.., CXorf..) are names of uncharacterised human loci and say nothing
# about an ortholog in another organism. NCBI also flags symbols containing 'like' (PC3LIKE). Both are
# NCBI BAD_GENE_NAME in a non-human submission.
BAD_SYMBOL = re.compile(r'^C(\d+|X|Y)orf\d+(_\d+)?$|like', re.I)

span_count = collections.Counter()
for b in blocks:
    if b[0] == 'feat' and b[1][0].split('\t')[2] == 'gene':
        span_count[(b[2], b[1][0].split('\t')[0].strip('<>'), b[1][0].split('\t')[1].strip('<>'))] += 1
cur_xref = None
for b in blocks:
    if b[0] != 'feat':
        out.extend(b[1])
        continue
    if b[1][0].split('\t')[2] == 'gene':
        g = b[1][0].split('\t')
        q = {l.split('\t')[3]: l.split('\t')[4] for l in b[1] if l.startswith('\t\t\t') and len(l.split('\t')) > 4}
        cur_xref = None
        if span_count[(b[2], g[0].strip('<>'), g[1].strip('<>'))] > 1:
            cur_xref = [('locus_tag', q.get('locus_tag'))] + \
                ([('gene', q['gene'])] if 'gene' in q and not BAD_SYMBOL.search(q['gene']) else [])
    ls = b[1]
    first = ls[0].split('\t')
    ftype = first[2]
    quals = [l for l in ls if l.startswith('\t\t\t')]
    product = next((l.split('\t')[4] for l in quals if l.split('\t')[3] == 'product'), None)
    new = [ls[0]]
    for l in ls[1:]:
        if l.startswith('\t\t\t'):
            q = l.split('\t')
            key, val = q[3], q[4] if len(q) > 4 else ''
            if key == 'gene' and BAD_SYMBOL.search(val):
                cnt[f'bad_symbol_removed_{ftype}'] += 1
                continue
            if key == 'product' and ftype == 'CDS' and re.search(r', variant \d+$', val):
                # isoforms of a gene share the protein name; 'variant' stays on the mRNA (GENE_PRODUCT_CONFLICT)
                val = re.sub(r', variant \d+$', '', val)
                cnt['CDS_product_variant_removed'] += 1
                l = f'\t\t\t{key}\t{val}'
            if key == 'EC_number' and val in EC_CHANGE:
                if EC_CHANGE[val] is None:
                    cnt['EC_ecn_removed'] += 1
                    continue
                val = EC_CHANGE[val]
                cnt['EC_ecn_replaced'] += 1
                l = f'\t\t\t{key}\t{val}'
            if key == 'EC_number':
                if product and product.lower().startswith('hypothetical protein'):
                    cnt['EC_removed_hypothetical'] += 1
                    continue
                parts = val.split('.')
                if len(parts) < 4:
                    val = '.'.join(parts + ['-'] * (4 - len(parts)))
                    cnt['EC_padded'] += 1
                    l = f'\t\t\t{key}\t{val}'
            if key == 'product' and re.search(r'[.,;:]\s*$', val):
                val = re.sub(r'[.,;:\s]+$', '', val)
                cnt[f'product_trailing_char_removed_{ftype}'] += 1
                l = f'\t\t\t{key}\t{val}'
        new.append(l)
    if ftype in ('mRNA', 'CDS') and cur_xref:
        last_q = max(i for i, x in enumerate(new) if x.startswith('\t\t\t') or i == 0)
        for k, v in reversed(cur_xref):
            if not any(l.split('\t')[3:4] == [k] for l in quals):
                new.insert(last_q + 1, f'\t\t\t{k}\t{v}')
        cnt[f'gene_xref_added_{ftype}'] += 1
    if ftype == 'ncRNA' and not any(l.split('\t')[3:4] == ['ncRNA_class'] for l in quals):
        a, z = int(first[0].lstrip('<>')), int(first[1].lstrip('<>'))
        strand = '+' if a <= z else '-'
        key = (b[2], min(a, z), max(a, z), strand)
        c = ncls.get(key)
        if c:
            last_q = max(i for i, x in enumerate(new) if x.startswith('\t\t\t') or i == 0)
            new.insert(last_q + 1, f'\t\t\tncRNA_class\t{c}')
            cnt['ncRNA_class_added'] += 1
        else:
            cnt['ncRNA_class_not_found'] += 1
    out.extend(new)

with open(tbl_out, 'w') as o:
    o.write('\n'.join(out))
with open(rep, 'w') as o:
    for k in sorted(cnt):
        o.write(f'{k}\t{cnt[k]}\n')
print(dict(cnt))
if cnt['ncRNA_class_not_found']:
    sys.exit(f"ERROR: {cnt['ncRNA_class_not_found']} ncRNA features without a class in the refined GFF3")
