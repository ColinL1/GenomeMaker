#!/usr/bin/env python3
"""Step 0d: apply contamination-screen actions and NCBI sequence hygiene, and lift the GFF3.

Sequence changes, in this order:
  1. actions from an FCS action table (FCS-GX *.fcs_gx_report.txt and/or FCS-adaptor
     fcs_adaptor_report.txt; columns accession, length, action, range(s)):
       EXCLUDE -> drop the sequence; TRIM -> remove the range (must touch a sequence end);
       TRIM_END (screen_review.py 'trim') -> remove the range and the bases between it and the nearer
       sequence end (the range must lie within 25 bp of that end);
       FIX -> remove the range if it touches an N run (a contig end next to a gap; the gap stays
       100 N), otherwise replace it with N. REVIEW/INFO are listed, not applied.
     Extra ranges (e.g. accepted VecScreen matches) can be given with --mask-bed (0-based, masked to N).
  2. terminal N runs are removed (NCBI does not accept sequences that begin or end with N);
  3. N runs are upper-cased ('n' from soft-masking inside gaps), other bases keep their case.
GFF3 lifting: coordinates are shifted by the bases removed upstream of them; a gene (with all its
children) is dropped when one of its exon/CDS/UTR/codon/RNA features overlaps a removed or masked
range, or when it lies on an excluded sequence; gene/mRNA spans that only cross a removed range
(an intron) are shortened. ##sequence-region lines are rewritten. Also writes deletions.tsv
(seq, start0, end0 of every removed range) and lengths.tsv (seq, new length) for lifting other files.

usage: apply_screen_actions.py --fasta in.fa --gff in.gff3 [--gff more.gff3 ...] --out-fasta out.fa
          --out-gff-dir DIR --out-log actions.tsv [--report fcs_report.txt ...] [--mask-bed x.bed]
"""
import argparse
import collections
import os
import re
import sys


def read_fasta(path):
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


ap = argparse.ArgumentParser()
ap.add_argument('--fasta', required=True)
ap.add_argument('--gff', action='append', default=[])
ap.add_argument('--report', action='append', default=[])
ap.add_argument('--mask-bed')
ap.add_argument('--out-fasta', required=True)
ap.add_argument('--out-gff-dir', required=True)
ap.add_argument('--out-log', required=True)
a = ap.parse_args()

actions = collections.defaultdict(list)   # seq -> [(action, start0, end0)]
log = []
for rp in a.report:
    for line in open(rp):
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) >= 8 and f[1].isdigit() and f[2].isdigit():
            # FCS-GX: seq_id, start_pos, end_pos, seq_len, action, div, agg_cont_cov, top_tax_name
            seq, act = f[0], f[4].upper()
            ranges = [(int(f[1]) - 1, int(f[2]))]
            f = [seq, f[3], act, f'{f[1]}..{f[2]}']
        else:
            # FCS-adaptor: accession, length, action, range(s)
            seq, act = f[0], f[2].replace('ACTION_', '').upper()
            ranges = []
            for r in f[3].split(','):
                if '..' in r:
                    x, y = r.split('..')
                    ranges.append((int(x) - 1, int(y)))
        if act in ('REVIEW', 'INFO', 'REVIEW_RARE'):
            log.append(('not_applied', seq, act, f[3], os.path.basename(rp)))
            continue
        if act == 'EXCLUDE':
            actions[seq].append(('EXCLUDE', 0, 0))
        else:
            for x, y in ranges:
                actions[seq].append((act, x, y))
        log.append(('applied', seq, act, f[3], os.path.basename(rp)))
if a.mask_bed:
    for line in open(a.mask_bed):
        f = line.rstrip('\n').split('\t')
        actions[f[0]].append(('FIX', int(f[1]), int(f[2])))
        log.append(('applied', f[0], 'FIX', f'{int(f[1]) + 1}..{f[2]}', 'mask_bed'))

newlen, excluded = {}, set()
deleted = collections.defaultdict(list)   # seq -> [(start0, end0)] removed from the sequence
masked = collections.defaultdict(list)    # seq -> [(start0, end0)] replaced by N
order = []
with open(a.out_fasta, 'w') as fo:
    for name, s in read_fasta(a.fasta):
        L = len(s)
        acts = actions.get(name, [])
        if any(x[0] == 'EXCLUDE' for x in acts):
            excluded.add(name)
            continue
        dels = []
        seq = list(s)
        for act, x, y in acts:
            if act == 'FIX':
                touches_gap = (x > 0 and s[x - 1] in 'Nn') or (y < L and s[y] in 'Nn')
                if touches_gap:
                    dels.append((x, y))
                    log.append(('applied', name, 'FIX_DELETED_AT_GAP', f'{x + 1}..{y}', 'next to an N run'))
                else:
                    for i in range(x, y):
                        seq[i] = 'N'
                    masked[name].append((x, y))
            elif act == 'TRIM':
                if x != 0 and y != L:
                    sys.exit(f'ERROR: TRIM range {name}:{x + 1}..{y} touches no sequence end')
                dels.append((x, y))
            elif act == 'TRIM_END':
                x2, y2 = (0, y) if x <= L - y else (x, L)
                if min(x, L - y) > 25:
                    sys.exit(f'ERROR: TRIM_END range {name}:{x + 1}..{y} is more than 25 bp from a sequence end')
                dels.append((x2, y2))
        s2 = ''.join(seq)
        # kept segments of the old sequence (0-based half-open), after the deletions
        segs, pos = [], 0
        for x, y in sorted(dels):
            if x > pos:
                segs.append((pos, x))
            pos = max(pos, y)
        if pos < L:
            segs.append((pos, L))
        body = ''.join(s2[x:y] for x, y in segs)

        def old_pos(k):
            """0-based index in body -> 0-based old position."""
            for x, y in segs:
                if k < y - x:
                    return x + k
                k -= y - x
            raise IndexError(k)
        m1 = len(body) - len(body.lstrip('Nn'))
        m2 = len(body) - len(body.rstrip('Nn'))
        if m1:
            dels.append((old_pos(0), old_pos(m1 - 1) + 1))
            log.append(('applied', name, 'TRIM_TERMINAL_N', f'{old_pos(0) + 1}..{old_pos(m1 - 1) + 1}', 'ncbi_hygiene'))
        if m2 and m2 < len(body):
            dels.append((old_pos(len(body) - m2), old_pos(len(body) - 1) + 1))
            log.append(('applied', name, 'TRIM_TERMINAL_N',
                        f'{old_pos(len(body) - m2) + 1}..{old_pos(len(body) - 1) + 1}', 'ncbi_hygiene'))
        core = body[m1:len(body) - m2]
        out = re.sub(r'[Nn]+', lambda m: m.group().upper(), core)
        # merge deletion intervals
        merged = []
        for x, y in sorted(dels):
            if merged and x <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], y)
            else:
                merged.append([x, y])
        deleted[name] = [tuple(m) for m in merged]
        newlen[name] = len(out)
        order.append(name)
        fo.write(f'>{name}\n')
        for i in range(0, len(out), 80):
            fo.write(out[i:i + 80] + '\n')


def lift(seq, pos):
    """1-based old position -> 1-based new position, or None if inside a deleted range."""
    shift = 0
    for x, y in deleted.get(seq, []):
        if pos - 1 >= y:
            shift += y - x
        elif pos - 1 >= x:
            return None
    return pos - shift


with open(os.path.join(os.path.dirname(os.path.abspath(a.out_fasta)), 'deletions.tsv'), 'w') as o:
    for n in order:
        for x, y in deleted[n]:
            o.write(f'{n}\t{x}\t{y}\n')
with open(os.path.join(os.path.dirname(os.path.abspath(a.out_fasta)), 'lengths.tsv'), 'w') as o:
    for n in order:
        o.write(f'{n}\t{newlen[n]}\n')

os.makedirs(a.out_gff_dir, exist_ok=True)
for gp in a.gff:
    recs = []
    parent = {}
    lines = [l.rstrip('\n') for l in open(gp)]
    top_of = {}
    feats = []
    for l in lines:
        if l.startswith('#') or not l.strip():
            continue
        f = l.split('\t')
        at = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
        feats.append((f, at))
        if 'ID' in at:
            parent[at['ID']] = at.get('Parent')

    def top(i):
        seen = 0
        while parent.get(i) and seen < 20:
            i = parent[i]
            seen += 1
        return i

    drop_top = set()
    for f, at in feats:
        key = at.get('ID') or at.get('Parent')
        t = top(key) if key else None
        if f[0] in excluded:
            drop_top.add(t)
            continue
        x, y = int(f[3]) - 1, int(f[4])
        span = f[2] in ('gene', 'mRNA', 'pseudogene', 'repeat_region')
        bad = [r for r in deleted.get(f[0], []) + masked.get(f[0], []) if x < r[1] and y > r[0]]
        if bad and not span:
            drop_top.add(t)
        elif span and (lift(f[0], int(f[3])) is None or lift(f[0], int(f[4])) is None):
            drop_top.add(t)
        elif span and any(x < r[1] and y > r[0] for r in masked.get(f[0], [])):
            drop_top.add(t)
    with open(os.path.join(a.out_gff_dir, os.path.basename(gp)), 'w') as fo:
        fo.write('##gff-version 3\n')
        used = {f[0] for f, _ in feats}
        for n in order:
            if n in used or 'sequence-region' in ''.join(lines[:3]):
                fo.write(f'##sequence-region {n} 1 {newlen[n]}\n')
        for f, at in feats:
            key = at.get('ID') or at.get('Parent')
            if key and top(key) in drop_top:
                continue
            f = list(f)
            f[3], f[4] = str(lift(f[0], int(f[3]))), str(lift(f[0], int(f[4])))
            fo.write('\t'.join(f) + '\n')
    for t in sorted(x for x in drop_top if x):
        log.append(('gff_dropped', t, 'overlaps removed/masked range or excluded sequence', '', os.path.basename(gp)))

with open(a.out_log, 'w') as o:
    o.write('status\tseq_or_feature\taction\trange\tsource\n')
    for r in log:
        o.write('\t'.join(r) + '\n')
c = collections.Counter((r[0], r[2]) for r in log)
for k in sorted(c):
    print(f'{k[0]}\t{k[1]}\t{c[k]}')
print(f'sequences in {len(order) + len(excluded)}, excluded {len(excluded)}, written {len(order)}')
