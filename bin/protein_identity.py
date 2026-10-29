#!/usr/bin/env python3
"""Step 3b: proteins of the refined GFF must equal the current set, except the transcripts in splits.tsv
(split pieces must be substrings of their source protein; status dropped: removed by refine_annotation;
status cds_edited: CDS corrected by an owner edit, not compared)
and partial proteins that refine_annotation extended by one residue at a partial end
(cds_3prime_restored, cds_5prime_to_seq_end).
usage: protein_identity.py refined.faa current.faa splits.tsv out.tsv [screen_actions.tsv]
Exits non-zero on any unexpected difference. Stop symbols ('*', '.') at the end are ignored.
With screen_actions.tsv (step 0d), transcripts of genes dropped by the screening are expected to be
missing and are counted separately."""
import sys
def fa(p):
    d, n, b = {}, None, []
    for l in open(p):
        if l.startswith('>'):
            if n: d[n] = ''.join(b)
            n, b = l[1:].split()[0], []
        else: b.append(l.strip())
    if n: d[n] = ''.join(b)
    return {k: v.upper().rstrip('*.') for k, v in d.items()}
new, cur = fa(sys.argv[1]), fa(sys.argv[2])
dropped_genes = set()
if len(sys.argv) > 5:
    for l in open(sys.argv[5]):
        f = l.rstrip('\n').split('\t')
        if f[0] == 'gff_dropped':
            dropped_genes.add(f[1])
gene_of = lambda t: t.rsplit('.', 1)[0]
screened_out = 0
old2new, dropped = {}, set()
for l in open(sys.argv[3]):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'gene': continue
    if f[4] == 'dropped': dropped.add(f[1])
    else: old2new.setdefault(f[1], []).append((f[3], f[4]))
split_old = set(old2new)
split_new = {n for v in old2new.values() for n, _ in v}
identical = changed_expected = restored_3prime = extended_5prime = 0; bad = []
for t, s in cur.items():
    if t in split_old or t in dropped:
        continue
    if t not in new:
        if gene_of(t) in dropped_genes:
            screened_out += 1
            continue
        bad.append((t, 'missing')); continue
    if new[t] != s:
        # refine_annotation cds_3prime_restored / cds_5prime_to_seq_end: a partial end extended by
        # 3 bp / 1 bp, so the protein gains at most one residue there (incomplete codons read as X)
        a_, b_ = s.strip('X'), new[t].strip('X')
        k = b_.find(a_)
        if k >= 0 and len(b_) > len(a_) and k <= 2 and len(b_) - len(a_) - k <= 2:
            restored_3prime += len(b_) - len(a_) - k > 0
            extended_5prime += k > 0
            continue
        bad.append((t, 'sequence differs')); continue
    identical += 1
moved_ok = edited = 0
for t, v in old2new.items():
    for n, st in v:
        if n not in new:
            if gene_of(n) in dropped_genes:
                screened_out += 1
                continue
            bad.append((n, 'split product missing')); continue
        if st == 'cds_edited':      # owner CDS correction (refine_annotation --transcript-edits)
            if n not in new: bad.append((n, 'edited transcript missing'))
            else: edited += 1
            continue
        if st == 'moved_unchanged':
            if new[n] != cur[t]: bad.append((n, f'moved_unchanged but differs from {t}'))
            else: moved_ok += 1
        else:
            changed_expected += 1
            if new[n] not in cur[t]: bad.append((n, f'split piece is not a substring of {t}'))
extra = set(new) - set(cur) - split_new
for t in extra: bad.append((t, 'unexpected new transcript'))
with open(sys.argv[4], 'w') as o:
    o.write(f'metric\tvalue\ncurrent_proteins\t{len(cur)}\nrefined_proteins\t{len(new)}\n'
            f'identical_unsplit\t{identical}\ncds_3prime_restored\t{restored_3prime}\n'
            f'cds_5prime_extended\t{extended_5prime}\n'
            f'split_source_transcripts\t{len(split_old)}\n'
            f'split_moved_unchanged\t{moved_ok}\nsplit_new_proteins\t{changed_expected}\n'
            f'cds_edited_by_owner\t{edited}\n'
            f'dropped_by_refine\t{len(dropped - split_old)}\n'
            f'removed_by_screening\t{screened_out}\nproblems\t{len(bad)}\n')
    for t, w in bad: o.write(f'PROBLEM\t{t}\t{w}\n')
print(open(sys.argv[4]).read())
sys.exit(1 if bad else 0)
