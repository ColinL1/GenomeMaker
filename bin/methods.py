#!/usr/bin/env python3
"""Methods text of the run, from what ran: tool versions and the automatic decisions.

The template (built by main.nf from the steps and settings the run used) holds two kinds of
placeholders:
  {{tool}}     the version(s) the run's processes reported for that versions.yml key
               ('a / b' when processes disagree; '(version not recorded)' when none did)
  {{=result}}  a sentence read from the run's decision files (assembly, polish, purge, mito,
               length, screen); '(not recorded)' when the file is absent
Nothing about results is written into the template itself.

usage: methods.py template.md all_versions.yml methods.md [--ranking assembly_ranking.tsv]
         [--polish-report polish_qc_report.tsv --polish-decision polish_decision.txt]
         [--purge-report purge_qc_report.tsv] [--mito-removed mito_contigs_removed.tsv]
         [--length-removed length_removed.tsv] [--screen-actions screen_actions.tsv]
Optional files may be NO_* placeholders.
"""
import argparse
import collections
import csv
import os
import re

ap = argparse.ArgumentParser()
ap.add_argument('template')
ap.add_argument('versions')
ap.add_argument('out')
for k in ('ranking', 'polish-report', 'polish-decision', 'purge-report', 'mito-removed',
          'length-removed', 'screen-actions'):
    ap.add_argument('--' + k)
a = ap.parse_args()


def usable(p):
    return bool(p) and os.path.exists(p) and not os.path.basename(p).startswith('NO_') and os.path.getsize(p) > 0


def rows(p):
    return list(csv.DictReader(open(p), delimiter='\t')) if usable(p) else []


vers = collections.defaultdict(set)
for l in open(a.versions, errors='replace'):
    m = re.match(r'^\s+([^:\s]+):\s*"?([^"]*?)"?\s*$', l.rstrip('\n'))
    if m and m.group(1) != 'container' and m.group(2) and m.group(2).lower() not in ('unknown', 'stub'):
        vers[m.group(1)].add(m.group(2))


def assembly():
    sel = [r for r in rows(a.ranking) if r.get('selected') == 'yes']
    if not sel:
        return None
    r = sel[0]
    dev = r.get('size_dev', 'NA')
    dev_txt = f', {float(dev):+.1%} from the expected size' if dev not in ('', 'NA') else ''
    return (f"{r['assembler']} selected by {r.get('decided_by', 'the ranking')} (BUSCO C {r['complete_pct']}%, "
            f"D {r['duplicated_pct']}%, {int(float(r['total_bp'])):,} bp{dev_txt})")


def polish():
    rep = {r['metric']: r for r in rows(a.polish_report)}
    if not rep:
        return 'the polished contigs were used without a regression check (run_polish_qc=false)'
    dec = open(a.polish_decision).read().strip() if usable(a.polish_decision) else '?'
    parts = []
    if 'merqury_qv' in rep:
        r = rep['merqury_qv']
        parts.append(f"Merqury QV {float(r['pre_polish']):.2f} -> {float(r['post_polish']):.2f}")
    if 'busco_completeness' in rep:
        r = rep['busco_completeness']
        parts.append(f"BUSCO C {float(r['pre_polish']):.1f}% -> {float(r['post_polish']):.1f}%")
    return f"the polish QC gate kept the {dec} contigs ({'; '.join(parts)})"


def purge():
    rep = {r['metric']: r for r in rows(a.purge_report)}
    if not rep:
        return None
    c, d = rep.get('busco_complete'), rep.get('busco_duplicated')
    out = []
    if c:
        out.append(f"BUSCO C {float(c['before_purge']):.1f}% -> {float(c['after_purge']):.1f}%")
    if d:
        out.append(f"D {float(d['before_purge']):.1f}% -> {float(d['after_purge']):.1f}%")
    verdicts = [r['verdict'] for r in rep.values()]
    gate = 'FAIL' if 'FAIL' in verdicts else ('WARN' if 'WARN' in verdicts else 'PASS')
    return f"{', '.join(out)} (purge QC gate {gate})"


def mito():
    if not a.mito_removed:
        return None
    rem = [r for r in rows(a.mito_removed) if r.get('action') == 'removed']
    if not rem:
        return 'no contig of the selected assembly matched the reference mitogenome'
    bp = sum(int(r['length']) for r in rem if r['length'].isdigit())
    return f"{len(rem)} contig(s) removed from the nuclear assembly ({', '.join(r['contig'] for r in rem)}; {bp:,} bp)"


def length():
    if not usable(a.length_removed):
        return None
    rem = rows(a.length_removed)
    return f"{len(rem)} sequence(s) removed ({sum(int(r['length']) for r in rem):,} bp)"


def screen():
    st = rows(a.screen_actions)
    if not usable(a.screen_actions):
        return None
    applied = collections.Counter(r['action'] for r in st if r.get('status') == 'applied')
    if not applied:
        return 'no action was needed'
    return 'actions applied: ' + ', '.join(f'{k} {v}' for k, v in sorted(applied.items()))


RESULTS = {'assembly': assembly, 'polish': polish, 'purge': purge, 'mito': mito, 'length': length, 'screen': screen}


def fill(m):
    key = m.group(1)
    if key.startswith('='):
        fn = RESULTS.get(key[1:])
        v = fn() if fn else None
        return v if v else '(not recorded)'
    v = vers.get(key)
    return ' / '.join(sorted(v)) if v else '(version not recorded)'


txt = re.sub(r'\{\{([^}]+)\}\}', fill, open(a.template).read())
open(a.out, 'w').write(txt)
print(txt)
