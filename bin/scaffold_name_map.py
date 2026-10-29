#!/usr/bin/env python3
"""WP10 / step 8: map final scaffold names to the input assembly names by sequence MD5.

Sequences are compared upper-cased, so soft-masking does not matter; a reverse-complement match is
reported as such. If the screening step (0d) changed sequences, pass the scaffolds before 0d with
--untrimmed (same names) and they are matched through that sequence.

usage: scaffold_name_map.py input_assembly.fa final_scaffolds.fa out.tsv [--untrimmed untrimmed.fa]
out.tsv: final_name, input_name, length, status (matched | reverse_complement | trimmed_substring | not_found),
         input_start (1-based, for trimmed_substring), strand
A final scaffold whose sequence (terminal N runs stripped) occurs inside one input sequence is
'trimmed_substring': purge_dups get_seqs -e removed contigs at its ends.
Input sequences with no final counterpart are listed with final_name '-' and status absent_in_final.
"""
import argparse
import hashlib

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


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


def md5(s):
    return hashlib.md5(s.upper().encode()).hexdigest()


ap = argparse.ArgumentParser()
ap.add_argument('input_fa')
ap.add_argument('final_fa')
ap.add_argument('out')
ap.add_argument('--untrimmed', help='final scaffolds before 0d trimming, same names')
a = ap.parse_args()

inp = {}
inp_len = {}
for n, s in read_fasta(a.input_fa):
    inp.setdefault(md5(s), n)
    inp_len[n] = len(s)
untr = dict(read_fasta(a.untrimmed)) if a.untrimmed else {}
found = set()
counts = {}
inp_seq = None
with open(a.out, 'w') as o:
    o.write('final_name\tinput_name\tlength\tstatus\tinput_start\tstrand\n')
    for n, s in read_fasta(a.final_fa):
        q = untr.get(n, s)
        h = md5(q)
        start, strand = '', '+'
        if h in inp:
            st, src = 'matched', inp[h]
        else:
            h2 = md5(q.translate(COMP)[::-1])
            if h2 in inp:
                st, src, strand = 'reverse_complement', inp[h2], '-'
            else:
                if inp_seq is None:
                    inp_seq = {k: v.upper() for k, v in read_fasta(a.input_fa)}
                core = q.upper().strip('N')
                st, src = 'not_found', '-'
                for k, v in inp_seq.items():
                    i = v.find(core)
                    if i >= 0:
                        st, src, start = 'trimmed_substring', k, str(i + 1)
                        break
                    j = v.find(core.translate(COMP)[::-1])
                    if j >= 0:
                        st, src, start, strand = 'trimmed_substring', k, str(j + 1), '-'
                        break
        if src != '-':
            found.add(src)
        counts[st] = counts.get(st, 0) + 1
        o.write(f'{n}\t{src}\t{len(s)}\t{st}\t{start}\t{strand}\n')
    for n, L in inp_len.items():
        if n not in found:
            o.write(f'-\t{n}\t{L}\tabsent_in_final\t\t\n')
            counts['absent_in_final'] = counts.get('absent_in_final', 0) + 1
print(counts)
