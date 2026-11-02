#!/usr/bin/env python3
"""purge_dups QC gate: BUSCO of the host bin before and after purge_dups.

FAIL (exit 1) when BUSCO completeness drops by more than --busco-tolerance percentage points: the
purge removed real sequence and the cutoffs need to be set by hand. WARN when duplication does not
drop: the purge removed no haplotigs that BUSCO can see.

usage: purge_qc_gate.py --before before.short_summary.txt --after after.short_summary.txt
         --cutoffs cutoffs --cutoffs-source 'calcuts (automatic)' --busco-tolerance 1.0
out: purge_qc_report.tsv, purge_qc_gate_mqc.tsv (MultiQC), purge_gate_status.txt (PASS | WARN: ...)
"""
import argparse
import re
import sys

ap = argparse.ArgumentParser()
ap.add_argument('--before', required=True)
ap.add_argument('--after', required=True)
ap.add_argument('--cutoffs', required=True)
ap.add_argument('--cutoffs-source', default='calcuts (automatic)')
ap.add_argument('--busco-tolerance', type=float, default=1.0)
a = ap.parse_args()


def busco(path):
    m = re.search(r'C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\]', open(path).read())
    if not m:
        sys.exit(f'ERROR: purge QC gate: no BUSCO summary line in {path}')
    return float(m.group(1)), float(m.group(3))


pre_c, pre_d = busco(a.before)
post_c, post_d = busco(a.after)
cut = open(a.cutoffs).read().split()
cut_txt = (f'lower {cut[0]}, upper {cut[5]} (calcuts line: {" ".join(cut)})' if len(cut) >= 6
           else (' '.join(cut) or 'empty cutoffs file'))

c_ok = post_c - pre_c >= -a.busco_tolerance
d_ok = post_d < pre_d
rows = [('busco_complete', f'{pre_c:.2f}', f'{post_c:.2f}', f'{post_c - pre_c:+.2f}', 'PASS' if c_ok else 'FAIL'),
        ('busco_duplicated', f'{pre_d:.2f}', f'{post_d:.2f}', f'{post_d - pre_d:+.2f}', 'PASS' if d_ok else 'WARN')]
header = 'metric\tbefore_purge\tafter_purge\tdelta\tverdict\n'
with open('purge_qc_report.tsv', 'w') as o:
    o.write(header)
    for r in rows:
        o.write('\t'.join(r) + '\n')

status = 'PASS' if c_ok and d_ok else ('FAIL' if not c_ok else 'WARN')
desc = (f'purge_dups gate: {status}. Cutoffs from {a.cutoffs_source}: {cut_txt}. '
        f'BUSCO completeness may drop by at most {a.busco_tolerance} pp.')
with open('purge_qc_gate_mqc.tsv', 'w') as m:
    m.write("# id: 'purge_qc_gate'\n")
    m.write("# section_name: '03b purge_dups QC gate'\n")
    m.write(f'# description: "{desc}"\n')
    m.write("# plot_type: 'table'\n")
    m.write(header)
    for r in rows:
        m.write('\t'.join(r) + '\n')

for r in rows:
    print('\t'.join(r))

if not c_ok:
    sys.stderr.write(
        f'\nERROR: purge_dups removed BUSCO genes: completeness {pre_c:.2f}% -> {post_c:.2f}% '
        f'({post_c - pre_c:+.2f} pp; tolerance {a.busco_tolerance} pp).\n'
        f'Cutoffs from {a.cutoffs_source}: {cut_txt}.\n'
        'Look at the coverage histogram (03b_purge_dups/PB.stat), set --purge_dups_lower, '
        '--purge_dups_transition and --purge_dups_upper, and re-run with -resume '
        '(or --run_purge_dups false).\n')
    sys.exit(1)

msg = 'PASS' if status == 'PASS' else (
    f'WARN: BUSCO duplication did not drop ({pre_d:.2f}% -> {post_d:.2f}%); cutoffs from {a.cutoffs_source}: {cut_txt}')
open('purge_gate_status.txt', 'w').write(msg + '\n')
print(f'\npurge_dups gate: {msg}')
