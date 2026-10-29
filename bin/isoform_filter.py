#!/usr/bin/env python3
"""Isoform evidence filter on the refined GFF3 (after refine_annotation.py, before protein extraction).

An alternative isoform is kept only with evidence of its own; every gene keeps at least one isoform.

  prep     proteins of every transcript of multi-isoform genes (for DIAMOND against the related
           proteomes) and the coverage regions (BED) of the end comparisons below.
  filter   the rules; writes the filtered GFF3, a TSV of every removed transcript with its rule and a
           splits.tsv with the removed transcripts as drops (a removed split piece: its own row).

End comparison: a transcript m whose span falls short of its gene's span at one end is compared with the
transcript c of the gene that reaches that end. Evidence: RNA-seq introns used by only one of the two
(BRAKER src=E intron hints, >= --min-reads), RNA-seq depth over c's exons beyond m's end (ext) against
m's own end exon (own), and homolog cover: the fraction of the best related protein covered by m vs c.
  m_supported        m has such an intron, or the homologs cover it better (>= --cov-diff)
  m_redundant        m's CDS lies inside c's (same frame), no own intron, homologs do not favour m
  c_supported_only   m adds sequence without evidence; c's end has an intron, reads (ext >= 0.3 own)
                     or homolog support
  c_end_doubtful     c's extension has no reads (ext < 0.05 own), no intron, homologs do not favour c
  undecided          anything else (low expression, conflicting signals)
For m_redundant pairs that differ only at the protein start, the related proteins vote on the start
(best HSP per species reaching an N-terminus, position within --start-tol aa): a short start supported by
>= 2 species and opposed by none keeps m.

Rules, applied once in this order:
  1. remove m when none of its ends is m_supported, c_end_doubtful or a supported short start
     (isoform_no_own_evidence, isoform_undecided);
  2. remove c of a c_end_doubtful end when c has no intron of its own with RNA-seq support
     (isoform_long_end_unsupported);
  3. short internal CDS exons (<= --short-exon bp) present in some but not all kept transcripts of a gene:
     the carriers are removed unless the exon is supported, by RNA-seq on both flanking introns or by
     the homologs (window of the exon residues +-13 aa: identity >= 0.8 and exon residues >= half
     identical in >= 2 species, none with window identity < 0.5) (isoform_short_exon_unsupported);
  4. a transcript partial (cut or truncated) at one end while a complete transcript of the gene reaches
     beyond that end gives a gene-partial conflict: the one with the lower homolog cover is removed, the
     complete one on a tie (isoform_partial_conflict);
  each rule keeps the transcript with the longest CDS when it would empty a gene.

usage: isoform_filter.py prep refined.gff3 genome.fa proteins.faa regions.bed
       isoform_filter.py filter refined.gff3 genome.fa hintsfile.gff.gz hits.tsv bedcov.tsv splits.tsv
                         out.gff3 removed.tsv out_splits.tsv summary.txt [options]
hits.tsv: DIAMOND -f 6 qseqid sseqid pident qstart qend sstart send qlen slen evalue bitscore qseq sseq,
subjects named <species>|<id>. bedcov.tsv: samtools bedcov -j regions.bed.
"""
import argparse
import collections
import gzip
import re
import sys

REVCOMP = str.maketrans('ACGTN', 'TGCAN')
CODE = dict(zip((a + b + c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'),
                'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))


def opener(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)


def attrs_of(s):
    return dict(kv.split('=', 1) for kv in s.strip().rstrip(';').split(';') if '=' in kv)


def read_fasta(p):
    seqs, name, buf = {}, None, []
    for line in open(p):
        if line.startswith('>'):
            if name is not None:
                seqs[name] = ''.join(buf).upper()
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(buf).upper()
    return seqs


def read_gff(p):
    """Genes with mRNAs: gene -> [transcript]; transcript -> model (chrom, strand, exons, cds, partial)."""
    genes, models = collections.OrderedDict(), {}
    for line in open(p):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        at = attrs_of(f[8])
        if f[2] == 'mRNA':
            genes.setdefault(at['Parent'], []).append(at['ID'])
            models[at['ID']] = dict(gene=at['Parent'], chrom=f[0], strand=f[6], start=int(f[3]), end=int(f[4]),
                                    exons=[], cds=[], partial=at.get('partial', ''))
        elif f[2] == 'exon' and at.get('Parent') in models:
            models[at['Parent']]['exons'].append((int(f[3]), int(f[4])))
        elif f[2] == 'CDS' and at.get('Parent') in models:
            models[at['Parent']]['cds'].append((int(f[3]), int(f[4]), int(f[7])))
    for m in models.values():
        m['exons'].sort()
        m['cds'].sort()
        m['introns'] = [(m['exons'][i][1] + 1, m['exons'][i + 1][0] - 1) for i in range(len(m['exons']) - 1)]
    return genes, models


def cds_positions(m):
    """Genomic positions of the CDS in transcript order, from the 5' phase on (whole codons)."""
    segs = sorted(m['cds'], reverse=(m['strand'] == '-'))
    pos = []
    for a, b, _ in segs:
        pos += list(range(a, b + 1)) if m['strand'] == '+' else list(range(b, a - 1, -1))
    return pos[segs[0][2]:] if segs else pos


def protein(m, seq):
    pos = cds_positions(m)
    nt = ''.join(seq[p - 1] for p in pos)
    if m['strand'] == '-':
        nt = nt.translate(REVCOMP)
    return ''.join(CODE.get(nt[i:i + 3], 'X') for i in range(0, len(nt) - 2, 3)).rstrip('*')


def frames(m):
    """{base: codon position}; equal values in two transcripts = same reading frame."""
    return {p: k % 3 for k, p in enumerate(cds_positions(m))}


def end_pairs(genes, models):
    """[(pair id, gene, m, c, '5p'/'3p', side 'L'/'R')] for every transcript end short of its gene's end;
    c = '-' when no transcript reaches beyond m on that side."""
    out, n = [], 0
    for g, ts in genes.items():
        if len(ts) < 2:
            continue
        lo, hi = min(models[t]['start'] for t in ts), max(models[t]['end'] for t in ts)
        for t in ts:
            m = models[t]
            for side in ('L', 'R'):
                if (side == 'L' and m['start'] <= lo) or (side == 'R' and m['end'] >= hi):
                    continue
                cands = [u for u in ts if u != t and (models[u]['start'] < m['start'] if side == 'L'
                                                      else models[u]['end'] > m['end'])]
                c = min(cands, key=lambda u: ((models[u]['start'] if side == 'L' else -models[u]['end']),
                                              -len(models[u]['exons']))) if cands else '-'
                end = ('5p' if m['strand'] == '+' else '3p') if side == 'L' else ('3p' if m['strand'] == '+' else '5p')
                n += 1
                out.append((f'p{n}', g, t, c, end, side))
    return out


def regions(pairs, models):
    """Coverage regions (BED, 0-based): pid|ext = c's exon bases beyond m's end; pid|own = m's end exon."""
    out = []
    for pid, g, t, c, end, side in pairs:
        if c == '-':
            continue
        m, mc = models[t], models[c]
        for a, b in mc['exons']:
            if side == 'L' and a < m['start']:
                out.append((m['chrom'], a - 1, min(b, m['start'] - 1), f'{pid}|ext'))
            if side == 'R' and b > m['end']:
                out.append((m['chrom'], max(a, m['end'] + 1) - 1, b, f'{pid}|ext'))
        a, b = m['exons'][0] if side == 'L' else m['exons'][-1]
        out.append((m['chrom'], a - 1, b, f'{pid}|own'))
    return out


def prep(a):
    genes, models = read_gff(a.gff)
    seqs = read_fasta(a.fasta)
    with open(a.faa, 'w') as o:
        for g, ts in genes.items():
            if len(ts) > 1:
                for t in ts:
                    o.write(f'>{t}\n{protein(models[t], seqs[models[t]["chrom"]])}\n')
    with open(a.bed, 'w') as o:
        for r in regions(end_pairs(genes, models), models):
            o.write('\t'.join(map(str, r)) + '\n')


# ----------------------------------------------------------------------------- filter
def filt(a):
    genes, models = read_gff(a.gff)
    seqs = read_fasta(a.fasta)
    reads = collections.Counter()
    for line in opener(a.hints):
        f = line.split('\t')
        if len(f) > 8 and f[2] == 'intron' and 'src=E' in f[8]:
            mm = re.search(r'mult=(\d+)', f[8])
            reads[(f[0], int(f[3]), int(f[4]))] += int(mm.group(1)) if mm else 1
    hsps = collections.defaultdict(list)          # query -> [(bits, species, subject, qs, ss, se, slen, qseq, sseq)]
    for line in open(a.hits):
        f = line.rstrip('\n').split('\t')
        if len(f) < 13:
            continue
        hsps[f[0]].append((float(f[10]), f[1].split('|')[0], f[1], int(f[3]), int(f[5]), int(f[6]), int(f[8]), f[11], f[12]))
    best_pair = {}                                 # (query, subject) -> best HSP
    for q, hs in hsps.items():
        for h in hs:
            k = (q, h[2])
            if k not in best_pair or h[0] > best_pair[k][0]:
                best_pair[k] = h
    cov = collections.defaultdict(lambda: [0.0, 0])
    for line in open(a.bedcov):
        f = line.rstrip('\n').split('\t')
        c = cov[f[3]]
        c[0] += float(f[4])
        c[1] += int(f[2]) - int(f[1])
    prot = {t: protein(m, seqs[m['chrom']]) for t, m in models.items() if len(genes[m['gene']]) > 1}

    def subject_cover(q, s):
        h = best_pair.get((q, s))
        return (h[5] - h[4] + 1) / h[6] if h else 0.0

    def best_subject(q):
        hs = hsps.get(q)
        return max(hs)[2] if hs else None

    def best_per_species(q):
        out = {}
        for h in sorted(hsps.get(q, []), reverse=True):
            out.setdefault(h[1], h)
        return out

    def own_supported(t, ts):
        m = models[t]
        others = {i for u in ts if u != t for i in models[u]['introns']}
        return [i for i in m['introns'] if i not in others and reads[(m['chrom'], *i)] >= a.min_reads]

    removed = {}                                   # transcript -> (rule, detail)
    kept = lambda g: [t for t in genes[g] if t not in removed]

    def guard(g, cands, rule, detail):
        """Remove cands from gene g, keeping its longest CDS if they would empty it."""
        rest = [t for t in kept(g) if t not in cands]
        cands = list(cands)
        if not rest and cands:
            keep = max(cands, key=lambda t: len(cds_positions(models[t])))
            cands.remove(keep)
        for t in cands:
            removed[t] = (rule, detail.get(t, '') if isinstance(detail, dict) else detail)

    # ---- end comparisons (rules 1 and 2)
    pairs = end_pairs(genes, models)
    per_m = collections.defaultdict(list)
    doubtful = []
    tally = collections.Counter()
    for pid, g, t, c, end, side in pairs:
        if c == '-':
            continue
        ts = genes[g]
        m_only = [i for i in models[t]['introns'] if i not in models[c]['introns']]
        c_only = [i for i in models[c]['introns'] if i not in models[t]['introns']]
        mr = max((reads[(models[t]['chrom'], *i)] for i in m_only), default=None)
        cr = max((reads[(models[c]['chrom'], *i)] for i in c_only), default=None)
        e, o = cov.get(pid + '|ext', [0, 0]), cov.get(pid + '|own', [0, 0])
        ed = e[0] / e[1] if e[1] else 0.0
        od = o[0] / o[1] if o[1] else 0.0
        m_j, c_j = mr is not None and mr >= a.min_reads, cr is not None and cr >= a.min_reads
        if m_j:
            rna = 'm_junction'
        elif c_j:
            rna = 'c_junction'
        elif od >= 2 and ed >= 0.3 * od:
            rna = 'reads_continue_into_c'
        elif od >= 2 and ed < 0.05 * od:
            rna = 'c_extension_no_reads'
        else:
            rna = 'low_or_unclear'
        sc, sm = best_subject(c), best_subject(t)
        if not sc and not sm:
            prot_v = 'no_homolog'
        else:
            d1 = subject_cover(c, sc) - subject_cover(t, sc) if sc else 0
            d2 = subject_cover(t, sm) - subject_cover(c, sm) if sm else 0
            prot_v = ('favours_c' if d1 >= a.cov_diff and d2 < a.cov_diff else
                      'favours_m' if d2 >= a.cov_diff and d1 < a.cov_diff else
                      'equal' if d1 < a.cov_diff and d2 < a.cov_diff else 'mixed')
        fm, fc = frames(models[t]), frames(models[c])
        contained = all(fc.get(x) == v for x, v in fm.items())
        if rna == 'm_junction' or prot_v == 'favours_m':
            v = 'm_supported'
        elif contained:
            v = 'm_redundant'
        elif rna in ('c_junction', 'reads_continue_into_c') or prot_v == 'favours_c':
            v = 'c_supported_only'
        elif rna == 'c_extension_no_reads' and prot_v != 'favours_c':
            v = 'c_end_doubtful'
        else:
            v = 'undecided'
        if v == 'm_redundant' and end == '5p' and prot[c].endswith(prot[t]):
            k = len(prot[c]) - len(prot[t])
            votes = collections.Counter()
            for sp, h in best_per_species(c).items():
                qs, ss = h[3], h[4]
                if qs > a.start_reach and ss > a.start_reach:
                    continue
                p = qs - ss + 1
                near_long, near_short = abs(p - 1) <= a.start_tol, abs(p - (k + 1)) <= a.start_tol
                if near_long and not near_short:
                    votes['long'] += 1
                elif near_short and not near_long:
                    votes['short'] += 1
            if votes['short'] >= 2 and not votes['long']:
                v = 'short_start_supported'
        tally[v] += 1
        per_m[(g, t)].append(v)
        if v == 'c_end_doubtful':
            doubtful.append((g, c))
    keep_v = {'m_supported', 'c_end_doubtful', 'short_start_supported'}
    by_gene = collections.defaultdict(dict)
    for (g, t), vs in per_m.items():
        if not keep_v & set(vs):
            by_gene[g][t] = ('isoform_undecided' if 'undecided' in vs else 'isoform_no_own_evidence', ','.join(sorted(set(vs))))
    for g, d in by_gene.items():
        for rule in ('isoform_no_own_evidence', 'isoform_undecided'):
            sel = [t for t, (r, _) in d.items() if r == rule]
            guard(g, sel, rule, {t: d[t][1] for t in sel})
    for g, c in doubtful:
        if c not in removed and not own_supported(c, genes[g]):
            guard(g, [c], 'isoform_long_end_unsupported', 'extension beyond another isoform has no reads, intron or homolog support')

    # ---- rule 3: short internal CDS exons
    def window_agrees(t, aa_from, aa_to):
        n = collections.Counter()
        for sp, hs in ((sp, [h for h in hsps.get(t, []) if h[1] == sp]) for sp in {h[1] for h in hsps.get(t, [])}):
            for h in sorted(hs, reverse=True):
                qpos, cols = h[3] - 1, []
                for x, y in zip(h[7], h[8]):
                    if x != '-':
                        qpos += 1
                    cols.append((qpos, x, y))
                qc = [p for p, x, _ in cols if x != '-']
                if not qc or min(qc) > aa_from - 13 or max(qc) < aa_to + 13:
                    continue
                w = [(x, y) for p, x, y in cols if aa_from - 13 <= p <= aa_to + 13]
                ex = [(x, y) for p, x, y in cols if aa_from <= p <= aa_to]
                wi = sum(x == y for x, y in w) / len(w)
                ei = sum(x == y for x, y in ex) / max(1, len(ex))
                n['agree' if wi >= 0.8 and ei >= 0.5 else 'disagree' if wi < 0.5 else 'partial'] += 1
                break
        return n['agree'] >= 2 and not n['disagree']

    for g in genes:
        ts = kept(g)
        if len(ts) < 2:
            continue
        short = collections.defaultdict(list)
        for t in ts:
            segs = sorted((x, y) for x, y, _ in models[t]['cds'])
            for x, y in segs[1:-1]:
                if y - x + 1 <= a.short_exon:
                    short[(x, y)].append(t)
        carriers, why = set(), {}
        for (x, y), car in short.items():
            if len(car) == len(ts):
                continue
            chrom = models[car[0]]['chrom']
            supported = False
            for t in car:
                ex = models[t]['exons']
                i = next((k for k, (ea, eb) in enumerate(ex) if ea <= x and y <= eb), None)
                if i is not None and 0 < i < len(ex) - 1:
                    left, right = (ex[i - 1][1] + 1, ex[i][0] - 1), (ex[i][1] + 1, ex[i + 1][0] - 1)
                    if ex[i] == (x, y) and reads[(chrom, *left)] >= a.min_reads and reads[(chrom, *right)] >= a.min_reads:
                        supported = True
                pos = cds_positions(models[t])
                idx = [k for k, p in enumerate(pos) if x <= p <= y]
                if idx and window_agrees(t, idx[0] // 3 + 1, idx[-1] // 3 + 1):
                    supported = True
            if not supported:
                for t in car:
                    carriers.add(t)
                    why[t] = f'short exon {x}-{y} unsupported; kept isoform without it exists'
        if carriers:
            guard(g, sorted(carriers), 'isoform_short_exon_unsupported', why)

    # ---- rule 4: partial end inside the span of a complete transcript
    def cover(t):
        s = best_subject(t)
        return subject_cover(t, s) if s else 0.0
    for g in genes:
        ts = kept(g)
        if len(ts) < 2:
            continue
        for t in ts:
            m = models[t]
            if t in removed or not m['partial']:
                continue
            p5 = m['partial'] in ('5prime', 'both')
            p3 = m['partial'] in ('3prime', 'both')
            left_p = p5 if m['strand'] == '+' else p3
            right_p = p3 if m['strand'] == '+' else p5
            for u in ts:
                if u == t or u in removed:
                    continue
                mu = models[u]
                u5 = mu['partial'] in ('5prime', 'both')
                u3 = mu['partial'] in ('3prime', 'both')
                u_left_p = u5 if mu['strand'] == '+' else u3
                u_right_p = u3 if mu['strand'] == '+' else u5
                clash = (left_p and not u_left_p and mu['start'] <= m['start']) or \
                        (right_p and not u_right_p and mu['end'] >= m['end'])
                if not clash:
                    continue
                loser = t if cover(t) <= cover(u) else u
                guard(g, [loser], 'isoform_partial_conflict',
                      f'partial {t} inside complete {u}; homolog cover {cover(t):.2f} vs {cover(u):.2f}')
                break

    # ---- outputs
    gene_span = {}
    for g in genes:
        ts = kept(g)
        if ts:
            gene_span[g] = (min(models[t]['start'] for t in ts), max(models[t]['end'] for t in ts))
    with open(a.gff) as fi, open(a.out_gff, 'w') as fo:
        for line in fi:
            if line.startswith('#'):
                fo.write(line)
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                fo.write(line)
                continue
            at = attrs_of(f[8])
            if f[2] == 'gene' and at.get('ID') in gene_span:
                f[3], f[4] = map(str, gene_span[at['ID']])
                fo.write('\t'.join(f) + '\n')
                continue
            ident = at.get('ID', '') if f[2] == 'mRNA' else at.get('Parent', '')
            if ident in removed:
                continue
            fo.write(line)
    with open(a.removed, 'w') as o:
        o.write('transcript\tgene\trule\tdetail\n')
        for t, (rule, det) in sorted(removed.items()):
            o.write(f'{t}\t{models[t]["gene"]}\t{rule}\t{det}\n')
    # a removed piece of a split transcript: its own row becomes a drop, so the source's other pieces
    # are still checked; any other removed transcript gets a row of its own
    with open(a.splits_in) as fi, open(a.out_splits, 'w') as o:
        pieces = set()
        for line in fi:
            f = line.rstrip('\n').split('\t')
            if len(f) > 4 and f[0] != 'gene' and f[3] in removed:
                pieces.add(f[3])
                line = '\t'.join([f[0], f[1], '', '', 'dropped', removed[f[3]][0]]) + '\n'
            o.write(line)
        for t, (rule, _) in sorted(removed.items()):
            if t not in pieces:
                o.write(f'{models[t]["gene"]}\t{t}\t\t\tdropped\t{rule}\n')
    with open(a.summary, 'w') as o:
        o.write(f'transcripts\t{len(models)}\nend_pairs\t{sum(1 for p in pairs if p[3] != "-")}\n')
        for k, n in sorted(tally.items()):
            o.write(f'end_verdict_{k}\t{n}\n')
        for k, n in sorted(collections.Counter(r for r, _ in removed.values()).items()):
            o.write(f'removed_{k}\t{n}\n')
        o.write(f'removed_total\t{len(removed)}\ntranscripts_kept\t{len(models) - len(removed)}\n')
    if any(not kept(g) for g in genes):
        sys.exit('ERROR: a gene lost all its transcripts')


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('prep')
    for k in ('gff', 'fasta', 'faa', 'bed'):
        p.add_argument(k)
    p = sub.add_parser('filter')
    for k in ('gff', 'fasta', 'hints', 'hits', 'bedcov', 'splits_in', 'out_gff', 'removed', 'out_splits', 'summary'):
        p.add_argument(k)
    p.add_argument('--min-reads', type=int, default=3)
    p.add_argument('--cov-diff', type=float, default=0.1)
    p.add_argument('--start-tol', type=int, default=5)
    p.add_argument('--start-reach', type=int, default=15)
    p.add_argument('--short-exon', type=int, default=15)
    a = ap.parse_args()
    prep(a) if a.cmd == 'prep' else filt(a)


if __name__ == '__main__':
    main()
