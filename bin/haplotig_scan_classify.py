#!/usr/bin/env python3
"""Haplotig scan, plus the regions purge_dups flagged.

  longest:  haplotig_scan_classify.py longest braker.gtf.gz braker.aa.gz longest.faa genes.tsv N_PLACED
            longest protein per BRAKER gene; genes.tsv = gene, transcript, chrom, start, end, placed, aa
  classify: haplotig_scan_classify.py classify genes.tsv allvsall.tsv out.tsv [flagged.bed]
            best near-identical copy on a placed scaffold (for placed genes: on another scaffold).
            flagged.bed (scaffold coordinates, 0-based, from purge_dups HAPLOTIG/OVLP)
            adds the fraction of each gene span inside flagged regions.

out.tsv: gene, chrom, placed, best_placed_copy, pident, qcov, class, purge_dups_flagged_frac
class: near_identical_placed_copy (pident >= 98, qcov >= 0.8), close_paralog_90, no_close_copy
"""
import collections
import gzip
import re
import sys


def opener(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)


def read_fasta(path):
    name, seq = None, []
    for line in opener(path):
        if line.startswith('>'):
            if name is not None:
                yield name, ''.join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name is not None:
        yield name, ''.join(seq)


def placed(chrom, n):
    m = re.match(r'scaffold_(\d+)$', chrom)
    return bool(m) and int(m.group(1)) <= n


def cmd_longest(gtf, aa, out_faa, out_tsv, n_placed):
    n_placed = int(n_placed)
    tx = {}
    for line in opener(gtf):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 9 and f[2] == 'transcript':
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            t = m.group(1) if m else f[8].strip()
            tx[t] = (t.rsplit('.', 1)[0], f[0], int(f[3]), int(f[4]))
    prot = {k: v.rstrip('*') for k, v in read_fasta(aa)}
    best = {}
    span = collections.defaultdict(lambda: [10 ** 12, 0])
    for t, (g, c, a, b) in tx.items():
        span[g][0] = min(span[g][0], a)
        span[g][1] = max(span[g][1], b)
        if t in prot and (g not in best or len(prot[t]) > len(prot[best[g]])):
            best[g] = t
    with open(out_faa, 'w') as fa, open(out_tsv, 'w') as o:
        o.write('gene\ttranscript\tchrom\tstart\tend\tplaced\taa\n')
        for g, t in best.items():
            c = tx[t][1]
            fa.write(f'>{g}\n{prot[t]}\n')
            o.write(f'{g}\t{t}\t{c}\t{span[g][0]}\t{span[g][1]}\t{int(placed(c, n_placed))}\t{len(prot[t])}\n')
    print(f'{len(best)} genes')


def cmd_classify(genes_tsv, allvsall, out_path, flagged=None):
    lab = {}
    for line in open(genes_tsv):
        f = line.rstrip('\n').split('\t')
        if f[0] != 'gene':
            lab[f[0]] = dict(chrom=f[2], start=int(f[3]), end=int(f[4]), placed=f[5] == '1')
    best = {}
    for line in open(allvsall):
        f = line.split('\t')
        q, s = f[0], f[1]
        if q == s or not lab[s]['placed']:
            continue
        if lab[q]['placed'] and lab[q]['chrom'] == lab[s]['chrom']:
            continue  # placed queries: ignore same-scaffold (tandem) copies
        pid, qcov = float(f[2]), (int(f[5]) - int(f[4]) + 1) / int(f[10])
        if q not in best or pid * qcov > best[q][0] * best[q][1]:
            best[q] = (pid, qcov, s)
    regions = collections.defaultdict(list)
    if flagged:
        for line in open(flagged):
            f = line.split('\t')
            regions[f[0]].append((int(f[1]) + 1, int(f[2])))
    cnt = collections.Counter()
    with open(out_path, 'w') as out:
        out.write('gene\tchrom\tplaced\tbest_placed_copy\tpident\tqcov\tclass\tpurge_dups_flagged_frac\n')
        for g, d in lab.items():
            b = best.get(g)
            if b and b[0] >= 98 and b[1] >= 0.8:
                c = 'near_identical_placed_copy'
            elif b and b[0] >= 90 and b[1] >= 0.5:
                c = 'close_paralog_90'
            else:
                c = 'no_close_copy'
            ov = sum(max(0, min(d['end'], e) - max(d['start'], s) + 1) for s, e in regions.get(d['chrom'], []))
            frac = ov / (d['end'] - d['start'] + 1)
            cnt[c] += 1
            cnt['purge_dups_flagged_ge50'] += frac >= 0.5
            out.write(f"{g}\t{d['chrom']}\t{int(d['placed'])}\t{b[2] if b else ''}\t{b[0] if b else ''}\t"
                      f"{round(b[1], 3) if b else ''}\t{c}\t{frac:.3f}\n")
    print(dict(cnt))


if __name__ == '__main__':
    {'longest': cmd_longest, 'classify': cmd_classify}[sys.argv[1]](*sys.argv[2:])
