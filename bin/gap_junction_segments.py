#!/usr/bin/env python3
"""Gap-fusion query builder: split BRAKER proteins at introns that contain a scaffold gap.

Every BRAKER transcript with an intron that fully contains an N run is cut into
protein segments at those introns. Each junction is labelled by intron-hint support:
'rna' (b2h/src=E hint at the exact intron), 'protein_only' (ProtHint only) or
'unsupported'. Supported junctions are the control set.

usage: gap_junction_segments.py --gtf braker.gtf.gz --aa braker.aa.gz \
          --hints hintsfile.gff.gz --fasta scaffolds.fa \
          --out-faa segments.faa --out-tsv junctions.tsv

junctions.tsv columns: transcript, junction, chrom, strand, intron_start, intron_end,
left_seg, right_seg, left_aa, right_aa, split_aa, support
"""
import argparse
import bisect
import collections
import gzip
import re

MIN_AA = 30


def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def read_fasta(path):
    name, seq = None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name is not None:
        yield name, ''.join(seq)


def braker_models(path):
    tx = {}
    with opener(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            if f[2] == 'transcript':
                m = re.search(r'transcript_id "([^"]+)"', f[8])
                t = m.group(1) if m else f[8].strip()
                tx[t] = dict(chrom=f[0], strand=f[6], cds=[], introns=[])
                continue
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            if not m:
                continue
            if f[2] == 'CDS':
                tx[m.group(1)]['cds'].append((int(f[3]), int(f[4])))
            elif f[2] == 'intron':
                tx[m.group(1)]['introns'].append((int(f[3]), int(f[4])))
    return tx


def intron_hints(path):
    h = collections.defaultdict(set)
    with opener(path) as fh:
        for line in fh:
            f = line.split('\t')
            if len(f) < 9 or f[2] != 'intron':
                continue
            s = re.search(r'src=([A-Z])', f[8])
            h[(f[0], int(f[3]), int(f[4]))].add(s.group(1) if s else '?')
    return h


def main():
    ap = argparse.ArgumentParser()
    for a in ('--gtf', '--aa', '--hints', '--fasta', '--out-faa', '--out-tsv'):
        ap.add_argument(a, required=True)
    a = ap.parse_args()

    gaps, starts = {}, {}
    for name, s in read_fasta(a.fasta):
        g = [(m.start() + 1, m.end()) for m in re.finditer(r'[Nn]+', s)]
        if g:
            gaps[name] = g
            starts[name] = [x[0] for x in g]

    def has_gap(c, lo, hi):
        if c not in gaps:
            return False
        i = bisect.bisect_left(starts[c], lo)
        return i < len(gaps[c]) and gaps[c][i][1] <= hi

    tx = braker_models(a.gtf)
    hints = intron_hints(a.hints)
    prot = {k: v.rstrip('*') for k, v in read_fasta(a.aa)}

    fa, tsv = open(a.out_faa, 'w'), open(a.out_tsv, 'w')
    tsv.write('transcript\tjunction\tchrom\tstrand\tintron_start\tintron_end\tleft_seg\tright_seg\t'
              'left_aa\tright_aa\tsplit_aa\tsupport\n')
    n_tx = n_j = 0
    for t, d in tx.items():
        gi = sorted(i for i in d['introns'] if has_gap(d['chrom'], *i))
        if not gi or t not in prot:
            continue
        p = prot[t]
        rev = d['strand'] == '-'
        cds = sorted(d['cds'], reverse=rev)
        splits = []
        for lo, hi in gi:
            off = sum(cb - ca + 1 for ca, cb in cds if (cb < lo if not rev else ca > hi))
            src = hints.get((d['chrom'], lo, hi), set())
            label = 'rna' if 'E' in src else ('protein_only' if src else 'unsupported')
            splits.append((off // 3, lo, hi, label))
        splits.sort()
        bounds = [0] + [s[0] for s in splits] + [len(p)]
        segs = []
        for k in range(len(bounds) - 1):
            sid = f'{t}|s{k}'
            s = p[bounds[k]:bounds[k + 1]]
            segs.append((sid, len(s)))
            if len(s) >= MIN_AA:
                fa.write(f'>{sid}\n{s}\n')
        n_tx += 1
        for k, (aa, lo, hi, label) in enumerate(splits):
            n_j += 1
            tsv.write(f'{t}\tj{k}\t{d["chrom"]}\t{d["strand"]}\t{lo}\t{hi}\t{segs[k][0]}\t{segs[k + 1][0]}\t'
                      f'{segs[k][1]}\t{segs[k + 1][1]}\t{aa}\t{label}\n')
    fa.close()
    tsv.close()
    print(f'{n_tx} transcripts with gap introns, {n_j} gap junctions')


if __name__ == '__main__':
    main()
