#!/usr/bin/env python3
"""Genome and annotation summary, computed from the delivered files only (final/ or release_curation/).

usage: genome_summary.py --fasta scaffolds.fa [--gff annotation.gff3 --stats stats.json --annotations annotations.txt]
         [--compleasm summary.txt] [--busco short_summary.txt] [--qv merqury.qv] [--repeat-tbl repeatmasker.tbl] [--hc-report release_sets_report.tsv]
         [--mito-fasta mitogenome.fasta] [--mito-check mito_gene_check.tsv] [--n-chrom 14]
         [--symbiont-busco short_summary.txt] [--funannotate-busco-db metazoa] --out summary.tsv
Lineage labels are read from the BUSCO summaries themselves; --funannotate-busco-db names the set
behind funannotate's 'busco' count. Without --gff/--stats (no annotation) only the assembly,
completeness, repeat and mitogenome rows are written.
out: section, metric, value
"""
import argparse
import collections
import json
import re

ap = argparse.ArgumentParser()
for k in ('fasta', 'gff', 'stats', 'annotations', 'compleasm', 'busco', 'qv', 'repeat-tbl', 'hc-report',
          'mito-fasta', 'mito-check', 'symbiont-busco', 'funannotate-busco-db', 'out'):
    ap.add_argument('--' + k)
ap.add_argument('--n-chrom', type=int, default=0, help='chromosome-scale scaffolds (scaffold_1..N); 0 = no row')
a = ap.parse_args()
rows = []


def add(section, metric, value):
    rows.append((section, metric, value))


def nx(lengths, frac=0.5):
    s = sorted(lengths, reverse=True)
    tot, acc = sum(s), 0
    for i, L in enumerate(s, 1):
        acc += L
        if acc >= tot * frac:
            return L, i
    return 0, 0


# ---------------- assembly
lens, contigs, gc = {}, [], collections.Counter()
lower = gaps = 0
name, buf = None, []


def seq_done():
    global lower, gaps
    if name is None:
        return
    s = ''.join(buf)
    lens[name] = len(s)
    gc.update(s.upper())
    lower += sum(1 for ch in s if ch.islower())
    pieces = re.split(r'[Nn]+', s)
    gaps += len(pieces) - 1
    contigs.extend(len(p) for p in pieces if p)


for line in open(a.fasta):
    if line.startswith('>'):
        seq_done()
        name, buf = line[1:].split()[0], []
    else:
        buf.append(line.strip())
seq_done()
total = sum(lens.values())
acgt = gc['A'] + gc['C'] + gc['G'] + gc['T']
s50, l50 = nx(lens.values())
c50, cl50 = nx(contigs)
chrom = [L for n, L in lens.items() if re.match(r'scaffold_(\d+)$', n) and int(n.split('_')[1]) <= a.n_chrom]
add('Assembly', 'Total length (bp)', f'{total:,}')
add('Assembly', 'Scaffolds', f'{len(lens):,}')
add('Assembly', 'Scaffold N50 (bp)', f'{s50:,}')
add('Assembly', 'Scaffold L50', l50)
add('Assembly', 'Largest scaffold (bp)', f'{max(lens.values()):,}')
if a.n_chrom:
    add('Assembly', f'Chromosome-scale scaffolds', f'{len(chrom)} ({sum(chrom):,} bp, {100 * sum(chrom) / total:.1f}%)')
add('Assembly', 'Contigs', f'{len(contigs):,}')
add('Assembly', 'Contig N50 (bp)', f'{c50:,}')
add('Assembly', 'Gaps', f'{gaps:,}')
add('Assembly', 'GC (%)', f'{100 * (gc["G"] + gc["C"]) / acgt:.2f}')
add('Assembly', 'N (%)', f'{100 * gc["N"] / total:.2f}')
add('Assembly', 'Soft-masked (%)', f'{100 * lower / total:.1f}')
if a.qv:
    f = open(a.qv).read().split()
    add('Assembly', 'Consensus QV (Merqury)', f'{float(f[3]):.2f}')
if a.compleasm:
    txt = open(a.compleasm).read()
    lin = re.search(r'lineage: (\S+)', txt)
    v = {k: (p, n) for k, p, n in re.findall(r'^([SDFIM]):([\d.]+)%, (\d+)', txt, re.M)}
    N = re.search(r'^N:(\d+)', txt, re.M)
    add('Completeness', f'compleasm {lin.group(1) if lin else ""} (n={N.group(1) if N else "?"})',
        f"C:{float(v['S'][0]) + float(v['D'][0]):.2f}% [S:{v['S'][0]}%, D:{v['D'][0]}%], F:{v['F'][0]}%, M:{v['M'][0]}%")


def busco_row(path, label):
    txt = open(path).read()
    lin = re.search(r'lineage dataset is: (\S+)', txt)
    m = re.search(r'C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)', txt.replace(' ', ''))
    if m:
        add('Completeness', f'{label} {lin.group(1) if lin else ""} (n={m.group(6)})',
            f'C:{m.group(1)}% [S:{m.group(2)}%, D:{m.group(3)}%], F:{m.group(4)}%, M:{m.group(5)}%')
    else:
        add('Completeness', label, txt.strip().splitlines()[0] if txt.strip() else 'NA')


if a.busco:
    busco_row(a.busco, 'BUSCO')
if a.symbiont_busco:
    busco_row(a.symbiont_busco, 'Symbiont bin BUSCO')

# ---------------- annotation (with an annotation)
if a.stats and a.gff:
    an = json.load(open(a.stats))['annotation']
    tl = an.get('transcript-level', {})
    pc = an['genes'] - an.get('tRNA', 0) - an.get('rRNA', 0) - an.get('ncRNA', 0)
    utr5, utr3, mrna_per_gene, gene_type = set(), set(), collections.Counter(), {}
    for line in open(a.gff):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        if f[2] in ('mRNA', 'tRNA', 'ncRNA', 'rRNA'):
            gene_type[re.search(r'Parent=([^;]+)', f[8]).group(1)] = f[2]
        if f[2] == 'mRNA':
            mrna_per_gene[re.search(r'Parent=([^;]+)', f[8]).group(1)] += 1
        elif f[2] in ('five_prime_UTR', 'three_prime_UTR'):
            (utr5 if f[2][0] == 'f' else utr3).add(re.search(r'Parent=([^;]+)', f[8]).group(1))
    add('Annotation', 'Protein-coding genes', pc)
    add('Annotation', 'mRNAs', an.get('mRNA'))
    add('Annotation', 'Genes with alternative isoforms', sum(1 for n in mrna_per_gene.values() if n > 1))
    ntype = collections.Counter(gene_type.values())     # RNA genes from the release GFF3 (rRNA fragments removed)
    add('Annotation', 'tRNA genes', ntype['tRNA'])
    add('Annotation', 'ncRNA genes', ntype['ncRNA'])
    add('Annotation', 'rRNA genes', ntype['rRNA'])
    add('Annotation', 'Mean gene length (bp)', round(an.get('avg_gene_length', 0)))
    add('Annotation', 'Mean protein length (aa)', round(tl.get('avg_protein_length', 0)))
    add('Annotation', 'Exons per mRNA (mean)', round(tl.get('total_cds_exons', 0) / max(1, an.get('mRNA', 1)), 1))
    add('Annotation', 'Single-exon mRNAs', tl.get('single_exon_transcript'))
    add('Annotation', 'mRNAs with complete CDS (start and stop codon)', tl.get('CDS_complete'))
    add('Annotation', 'mRNAs with partial CDS (no start, no stop, or neither)',
        tl.get('CDS_no-start', 0) + tl.get('CDS_no-stop', 0) + tl.get('CDS_no-start_no-stop', 0))
    add('Annotation', "mRNAs with 5' UTR", len(utr5))
    add('Annotation', "mRNAs with 3' UTR", len(utr3))
    if a.hc_report:
        d = dict(l.rstrip('\n').split('\t')[:2] for l in open(a.hc_report) if '\t' in l)
        add('Annotation', 'High-confidence protein-coding genes (subset)', d.get('hc_protein_coding_genes'))
    # product per gene from the release GFF3 (first mRNA); annotations.txt omits genes that share a span
    named = hyp = 0
    seen = set()
    for line in open(a.gff):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'mRNA':
            continue
        g = re.search(r'Parent=([^;]+)', f[8]).group(1)
        if g in seen:
            continue
        seen.add(g)
        m = re.search(r'product=([^;]+)', f[8])
        if m and not m.group(1).startswith('hypothetical protein'):
            named += 1
        else:
            hyp += 1
    add('Functional', 'Genes with a product name', named)
    add('Functional', 'Genes named hypothetical protein', hyp)
    add('Functional', 'Genes with a gene symbol', an.get('common_name'))
    labels = {'go_terms': 'mRNAs with GO terms', 'interproscan': 'mRNAs with InterPro domains', 'eggnog': 'mRNAs with eggNOG',
              'pfam': 'mRNAs with Pfam domains', 'cazyme': 'mRNAs with CAZyme', 'merops': 'mRNAs with MEROPS peptidase',
              'busco': f"mRNAs matching BUSCO ({a.funannotate_busco_db or 'funannotate'})",
              'secretion': 'mRNAs with secretion signal'}
    for k, v in tl.get('functional', {}).items():
        add('Functional', labels.get(k, k), v)

# ---------------- repeats
if a.repeat_tbl:
    # class rows start in column 0; their subfamily rows are indented; the class line itself counts
    # only elements without a subfamily, so class totals add the indented rows
    cls_bp, order, cur, nonN, tot_pct = collections.Counter(), [], None, None, None
    for line in open(a.repeat_tbl):
        if line.startswith('Repeat Stats'):
            break
        m = re.match(r'Total Length: (\d+)', line)
        if m:
            nonN = int(m.group(1))
        m = re.match(r'^(\s*)(\S+)\s+(\d+|--)\s+(\d+|--)\s+([\d.]+%|--)', line)
        if not m or m.group(2) in ('Class', 'ID', '====='):
            continue
        name_, bp = m.group(2), (int(m.group(4)) if m.group(4).isdigit() else 0)
        if name_ == 'Total':
            tot_pct = m.group(5)
            continue
        if not m.group(1):
            cur = name_
            order.append(cur)
        cls_bp[cur] += bp
    for c in sorted(order, key=lambda c: -cls_bp[c]):
        add('Repeats', f'{c} (% of non-N sequence)', f'{100 * cls_bp[c] / nonN:.2f}%')
    add('Repeats', 'Total masked (% of non-N sequence)', tot_pct)

# ---------------- mitogenome
if a.mito_fasta:
    L = sum(len(l.strip()) for l in open(a.mito_fasta) if not l.startswith('>'))
    add('Mitogenome', 'Length (bp)', f'{L:,}')
if a.mito_check:
    st = [l.rstrip('\n').split('\t') for l in open(a.mito_check)][1:]
    add('Mitogenome', 'Protein-coding genes full length', f"{sum(1 for r in st if r[-1] == 'OK')} of {len(st)}")

with open(a.out, 'w') as o:
    o.write('section\tmetric\tvalue\n')
    for r in rows:
        o.write('\t'.join(str(x) for x in r) + '\n')
print(open(a.out).read())
