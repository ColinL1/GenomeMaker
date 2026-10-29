#!/usr/bin/env python3
"""Step 9b (WP11, D3): relabel repeat classes after RepeatClassifier; the masked sequence is unchanged.

  1. Family labels come from the library headers ('>name#Class/Family'). Families labelled
     Unknown/unknown take the RepeatClassifier result (its *.classified FASTA).
  2. All labels are normalised (WP11.2): 'unknown' -> 'Unknown', 'X/unknown' -> 'X/Unknown',
     'DNAnona[/X]' and 'DNAauto[/X]' -> 'DNA[/X]'.
  3. Column 11 (class/family) of repeatmasker.out is rewritten from column 10 (family name).
     Simple_repeat, Low_complexity and other RepeatMasker-internal labels of families that are not
     in the library are kept as they are (after normalisation).

usage: repeat_relabel.py library.fa classified.fa repeatmasker.out out_library.fa out.out out_report.tsv
       repeat_relabel.py library.fa classified.fa - out_library.fa - out_report.tsv   (library only)
"""
import collections
import re
import sys

lib_p, cls_p, rmout_p, out_lib, out_out, out_rep = sys.argv[1:7]


def norm(label):
    lab = label
    if lab == 'unknown':
        lab = 'Unknown'
    lab = re.sub(r'/unknown$', '/Unknown', lab)
    lab = re.sub(r'^DNA(nona|auto)(?=/|$)', 'DNA', lab)
    return lab


def is_unknown(label):
    return label.split('/')[0].lower() == 'unknown'


lib = collections.OrderedDict()
cur = None
for line in open(lib_p):
    if line.startswith('>'):
        name, _, rest = line[1:].rstrip('\n').partition('#')
        label = rest.split()[0] if rest else 'Unknown'
        cur = name
        lib[cur] = [label, []]
    else:
        lib[cur][1].append(line)
new = {}
for line in open(cls_p):
    if line.startswith('>'):
        name, _, rest = line[1:].rstrip('\n').partition('#')
        new[name.split()[0]] = rest.split()[0] if rest else 'Unknown'

cnt = collections.Counter()
label = {}
for name, (lab, _) in lib.items():
    l2 = lab
    if is_unknown(lab):
        cnt['families_unknown_before'] += 1
        if name in new:
            l2 = new[name]
            cnt['families_classifier_result'] += 1
        else:
            cnt['families_missing_from_classifier'] += 1
    l3 = norm(l2)
    if l3 != l2:
        cnt['labels_normalised'] += 1
    if is_unknown(l3):
        cnt['families_unknown_after'] += 1
    label[name] = l3
with open(out_lib, 'w') as o:
    for name, (_, seq) in lib.items():
        o.write(f'>{name}#{label[name]}\n')
        o.writelines(seq)

bp_before, bp_after = collections.Counter(), collections.Counter()
tok = re.compile(r'\S+')
if rmout_p == '-':
    # library only (classification before RepeatMasker): no .out to relabel
    with open(out_rep, 'w') as o:
        o.write('metric\tvalue\n')
        for k in sorted(cnt):
            o.write(f'{k}\t{cnt[k]}\n')
    print(open(out_rep).read())
    sys.exit(0)
with open(rmout_p) as fi, open(out_out, 'w') as fo:
    for line in fi:
        toks = list(tok.finditer(line))
        if len(toks) < 11 or not toks[0].group().isdigit():
            fo.write(line)
            continue
        fam, old = toks[9].group(), toks[10].group()
        newlab = label.get(fam, norm(old))
        L = int(toks[6].group()) - int(toks[5].group()) + 1
        bp_before[old.split('/')[0]] += L
        bp_after[newlab.split('/')[0]] += L
        if newlab != old:
            cnt['out_records_relabelled'] += 1
            s, e = toks[10].span()
            pad = max(0, (e - s) - len(newlab))
            line = line[:s] + newlab + ' ' * pad + line[e:]
        cnt['out_records'] += 1
        fo.write(line)

tot = sum(bp_before.values())
with open(out_rep, 'w') as o:
    o.write('metric\tvalue\n')
    for k in sorted(cnt):
        o.write(f'{k}\t{cnt[k]}\n')
    for k in sorted(set(bp_before) | set(bp_after), key=lambda k: -bp_after.get(k, 0)):
        o.write(f'bp_{k}\t{bp_before.get(k, 0)}\t{bp_after.get(k, 0)}\n')
    ub = sum(v for k, v in bp_before.items() if k.lower() == 'unknown')
    ua = sum(v for k, v in bp_after.items() if k.lower() == 'unknown')
    o.write(f'unknown_share_of_repeat_bp_before\t{ub / tot:.4f}\nunknown_share_of_repeat_bp_after\t{ua / tot:.4f}\n')
print(open(out_rep).read())
