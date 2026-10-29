#!/usr/bin/env python3
"""Contamination screen review: list the hits that need a person, apply the person's decisions.

Clear-cut actions (FCS-GX EXCLUDE/TRIM/FIX, FCS-adaptor actions) are applied by
apply_screen_actions.py without review. The unclear ones need a decision:
  - FCS-GX REVIEW / REVIEW_RARE rows;
  - VecScreen Strong and Moderate matches (vecscreen_classify.py output).

decisions.tsv (optional; '#' lines and a 'seq' header are skipped):
  seq, start, end, decision, reason      decision = exclude | mask | trim | keep
  start/end are 1-based, as in the review list; 'exclude' drops the whole sequence, 'mask' replaces the
  range with N, 'trim' removes the range and the bases between it and the nearer sequence end (only
  within 25 bp of an end), 'keep' leaves it as it is.

Outputs
  screen_review.tsv       every item needing a decision, with the decision found (or 'UNDECIDED')
  review_excluded.tsv     FCS-GX-format rows (EXCLUDE, TRIM_END) for apply_screen_actions.py --report
  review_mask.bed         0-based ranges for apply_screen_actions.py --mask-bed
  accepted_vecscreen.tsv  kept VecScreen matches (seq, start, end, category, reason), for the QC gate
  accepted_fcs_gx.tsv     kept FCS-GX REVIEW rows (seq, start, end, action, reason), for the QC gate
Exit status 3 when any item is undecided, so the run stops for review.

usage: screen_review.py --vecscreen vecscreen_regions.tsv [--fcs-gx x.fcs_gx_report.txt ...]
         [--decisions decisions.tsv] --outdir DIR
"""
import argparse
import os
import sys

ap = argparse.ArgumentParser()
ap.add_argument('--fcs-gx', action='append', default=[])
ap.add_argument('--vecscreen', action='append', default=[])
ap.add_argument('--decisions')
ap.add_argument('--outdir', required=True)
a = ap.parse_args()


def real(p):
    return p and os.path.exists(p) and not os.path.basename(p).startswith('NO_')


items = []   # (seq, start, end, source, category, detail)
for p in filter(real, a.fcs_gx):
    for line in open(p):
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) >= 5 and f[4].upper() in ('REVIEW', 'REVIEW_RARE'):
            items.append((f[0], f[1], f[2], 'FCS-GX', f[4], ' '.join(f[5:8])))
for p in filter(real, a.vecscreen):
    for line in open(p):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'seq' or len(f) < 4:
            continue
        if f[3] in ('Strong', 'Moderate'):
            items.append((f[0], f[1], f[2], 'VecScreen', f[3], ' '.join(f[4:])))

dec = {}
if real(a.decisions):
    for line in open(a.decisions):
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        if f[0] == 'seq':
            continue
        d = f[3].strip().lower()
        if d not in ('exclude', 'mask', 'trim', 'keep'):
            sys.exit(f'ERROR: decisions: "{f[3]}" for {f[0]}:{f[1]}-{f[2]} is not exclude, mask, trim or keep')
        dec[(f[0], f[1], f[2])] = (d, f[4] if len(f) > 4 else '')

o = a.outdir
os.makedirs(o, exist_ok=True)
undecided = 0
with open(f'{o}/screen_review.tsv', 'w') as rv, open(f'{o}/review_excluded.tsv', 'w') as ex, \
        open(f'{o}/review_mask.bed', 'w') as mb, open(f'{o}/accepted_vecscreen.tsv', 'w') as av, \
        open(f'{o}/accepted_fcs_gx.tsv', 'w') as ag:
    rv.write('seq\tstart\tend\tsource\tcategory\tdetail\tdecision\treason\n')
    ex.write('#seq_id\tstart_pos\tend_pos\tseq_len\taction\tdiv\tagg_cont_cov\ttop_tax_name\n')
    av.write('seq\tstart\tend\tcategory\treason\n')
    ag.write('seq\tstart\tend\taction\treason\n')
    for seq, s, e, src, cat, det in items:
        d, why = dec.get((seq, s, e), ('UNDECIDED', ''))
        undecided += d == 'UNDECIDED'
        rv.write(f'{seq}\t{s}\t{e}\t{src}\t{cat}\t{det}\t{d}\t{why}\n')
        if d == 'exclude':
            ex.write(f'{seq}\t{s}\t{e}\t0\tEXCLUDE\treview\t0\t{why}\n')
        elif d == 'trim':
            ex.write(f'{seq}\t{s}\t{e}\t0\tTRIM_END\treview\t0\t{why}\n')
        elif d == 'mask':
            mb.write(f'{seq}\t{int(s) - 1}\t{e}\n')
        elif d == 'keep':
            (av if src == 'VecScreen' else ag).write(f'{seq}\t{s}\t{e}\t{cat}\t{why}\n')
print(f'{len(items)} items need a decision; undecided {undecided}')
if undecided:
    print(f'ERROR: {undecided} screen hits need a decision. Copy the UNDECIDED rows of screen_review.tsv '
          f'into the file given by --screen_decisions (columns seq, start, end, decision, reason; decision = '
          f'exclude, mask, trim or keep) and resume the run.', file=sys.stderr)
    sys.exit(3)
