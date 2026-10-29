#!/usr/bin/env python3
"""Drop sequences shorter than a minimum length (NCBI rejects sequences under 200 bp).

usage: length_filter.py in.fasta out.fasta removed.tsv MIN_LEN
removed.tsv: name, length
"""
import sys

inp, out_p, rem_p, min_len = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])


def records(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith('>'):
            if name is not None:
                yield name, ''.join(buf)
            name, buf = line[1:].rstrip('\n'), []
        else:
            buf.append(line.strip())
    if name is not None:
        yield name, ''.join(buf)


kept = removed = removed_bp = 0
with open(out_p, 'w') as o, open(rem_p, 'w') as r:
    r.write('name\tlength\n')
    for header, seq in records(inp):
        if len(seq) < min_len:
            r.write(f'{header.split()[0]}\t{len(seq)}\n')
            removed += 1
            removed_bp += len(seq)
            continue
        o.write(f'>{header}\n')
        for i in range(0, len(seq), 80):
            o.write(seq[i:i + 80] + '\n')
        kept += 1
print(f'kept {kept}; removed {removed} sequences shorter than {min_len} bp ({removed_bp} bp)')
