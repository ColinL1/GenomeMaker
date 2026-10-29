#!/usr/bin/env python3
"""Step 10 (WP5, D5): release GFF with tags, high-confidence subset, candidate lists.

funannotate does not carry custom GFF attributes, so the refined models are mapped onto the final
(funannotate) models by structure: (scaffold, strand, CDS segments, exon segments), falling back to
(scaffold, strand, CDS segments). Then:
  - the release GFF3 = final GFF3 with each mRNA's span set to its exon extent (funannotate writes
    the gene span on every isoform) + gene attributes evidence=, possible_haplotig_of=<final gene>,
    possible_haplotig=stageA and mRNA attribute note=possible_gap_fusion|possible_gap_misorder;
  - the high-confidence set: protein-coding genes with evidence != none, or with an InterPro, eggNOG
    or Pfam hit in annotations.txt; all tRNA/ncRNA/rRNA genes are kept;
  - lists: possible_haplotigs.tsv, gap_fusion_candidates.tsv, gap_misorder_candidates.tsv,
    splits.tsv (refined and final IDs), id_map.tsv.

usage: release_sets.py --final-gff F.gff3 --final-proteins F.proteins.fa --annotations F.annotations.txt
          --refined-gff refined.gff3 --tags gene_tags.tsv --splits splits.tsv --prefix Porites_sp --outdir DIR
"""
import argparse
import collections
import os
import sys


def attrs(s):
    return dict(kv.split('=', 1) for kv in s.strip().rstrip(';').split(';') if '=' in kv)


def models(path):
    """mRNA id -> dict(gene, chrom, strand, cds, exons); plus gene -> type of child."""
    mr, cds, ex = {}, collections.defaultdict(list), collections.defaultdict(list)
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        at = attrs(f[8])
        if f[2] == 'mRNA':
            mr[at['ID']] = dict(gene=at['Parent'], chrom=f[0], strand=f[6])
        elif f[2] == 'CDS':
            for p in at['Parent'].split(','):
                cds[p].append((int(f[3]), int(f[4])))
        elif f[2] == 'exon':
            for p in at['Parent'].split(','):
                ex[p].append((int(f[3]), int(f[4])))
    for t, d in mr.items():
        d['cds'] = tuple(sorted(cds[t]))
        d['exons'] = tuple(sorted(ex[t]))
    return mr


def main():
    ap = argparse.ArgumentParser()
    for k in ('--final-gff', '--final-proteins', '--annotations', '--refined-gff', '--tags', '--splits',
              '--prefix', '--outdir'):
        ap.add_argument(k, required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    fin = models(a.final_gff)
    ref = models(a.refined_gff)
    k_full = {(d['chrom'], d['strand'], d['cds'], d['exons']): t for t, d in fin.items()}
    k_cds = collections.defaultdict(list)
    for t, d in fin.items():
        k_cds[(d['chrom'], d['strand'], d['cds'])].append(t)
    tx_map, problems = {}, []
    for t, d in ref.items():
        f = k_full.get((d['chrom'], d['strand'], d['cds'], d['exons']))
        if not f:
            c = k_cds.get((d['chrom'], d['strand'], d['cds']), [])
            f = c[0] if len(c) == 1 else None
        if f:
            tx_map[t] = f
        else:
            problems.append(t)
    gene_map = {}
    for t, f in tx_map.items():
        gene_map.setdefault(ref[t]['gene'], fin[f]['gene'])
    with open(f'{a.outdir}/id_map.tsv', 'w') as o:
        o.write('refined_gene\trefined_transcript\tfinal_gene\tfinal_transcript\n')
        for t in sorted(ref):
            o.write(f"{ref[t]['gene']}\t{t}\t{gene_map.get(ref[t]['gene'], '')}\t{tx_map.get(t, '')}\n")

    tags = {}
    gene_tags = {}
    for line in open(a.tags):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'gene':
            continue
        g, t, ev, hof, hst, note = f
        tags[t] = note
        gene_tags[g] = (ev, hof, hst)

    # functional hits per final transcript (annotations.txt columns PFAM, InterPro, EggNog)
    hit_tx = set()
    hd = None
    for line in open(a.annotations):
        f = line.rstrip('\n').split('\t')
        if hd is None:
            hd = {k: i for i, k in enumerate(f)}
            continue
        if any(f[hd[c]].strip() for c in ('PFAM', 'InterPro', 'EggNog') if c in hd):
            hit_tx.add(f[hd['TranscriptID']])
    # annotations.txt writes transcripts as <id>-T<n>-style or the GFF ID; accept both
    fin_by_alt = {}
    for t in fin:
        fin_by_alt[t] = t
    hit_genes = {fin[fin_by_alt[t]]['gene'] for t in hit_tx if t in fin_by_alt}
    hit_genes |= {t.rsplit('-T', 1)[0] for t in hit_tx}

    final_gene_attr = collections.defaultdict(dict)
    final_tx_note = {}
    for rg, (ev, hof, hst) in gene_tags.items():
        fg = gene_map.get(rg)
        if not fg:
            continue
        at = final_gene_attr[fg]
        at['evidence'] = ev
        if hof:
            at['possible_haplotig_of'] = gene_map.get(hof, hof)
        if hst:
            at['possible_haplotig'] = hst
    for t, note in tags.items():
        if note and t in tx_map:
            final_tx_note[tx_map[t]] = note

    # exon extent per final mRNA (funannotate gives every isoform the gene's span)
    ext = {}
    for t, d in fin.items():
        if d['exons']:
            ext[t] = (d['exons'][0][0], max(y for _, y in d['exons']))
    n_span_fixed = 0

    hc_genes = set()
    gene_type = {}
    out_lines = []
    for line in open(a.final_gff):
        if line.startswith('#') or not line.strip():
            out_lines.append(line)
            continue
        f = line.rstrip('\n').split('\t')
        at = attrs(f[8])
        if f[2] == 'gene' and at['ID'] in final_gene_attr:
            extra = ';'.join(f'{k}={v}' for k, v in final_gene_attr[at['ID']].items())
            f[8] = f[8].rstrip(';') + ';' + extra
        elif f[2] == 'mRNA' and at['ID'] in final_tx_note:
            f[8] = f[8].rstrip(';') + ';note=' + final_tx_note[at['ID']]
        if f[2] == 'mRNA' and at['ID'] in ext and (int(f[3]), int(f[4])) != ext[at['ID']]:
            f[3], f[4] = str(ext[at['ID']][0]), str(ext[at['ID']][1])
            n_span_fixed += 1
        if f[2] in ('mRNA', 'tRNA', 'ncRNA', 'rRNA') and 'Parent' in at:
            gene_type[at['Parent']] = f[2]
        out_lines.append('\t'.join(f) + '\n')
    for g, typ in gene_type.items():
        if typ != 'mRNA':
            hc_genes.add(g)
        elif final_gene_attr.get(g, {}).get('evidence', 'none') != 'none' or g in hit_genes:
            hc_genes.add(g)
    with open(f'{a.outdir}/{a.prefix}.gff3', 'w') as o:
        o.writelines(out_lines)
    keep_ids = set()
    with open(f'{a.outdir}/{a.prefix}.hc.gff3', 'w') as o:
        for line in out_lines:
            if line.startswith('#') or not line.strip():
                o.write(line)
                continue
            f = line.rstrip('\n').split('\t')
            at = attrs(f[8])
            g = at['ID'] if f[2] == 'gene' else None
            if g is None:
                p = at.get('Parent', '')
                g = fin[p]['gene'] if p in fin else p
            if g in hc_genes:
                o.write(line)
                if f[2] == 'mRNA':
                    keep_ids.add(at['ID'])
    n = 0
    with open(a.final_proteins) as fi, open(f'{a.outdir}/{a.prefix}.hc.proteins.fa', 'w') as o:
        keep = False
        for line in fi:
            if line.startswith('>'):
                keep = line[1:].split()[0] in keep_ids
                n += keep
            if keep:
                o.write(line)

    def write_list(fn, rows, head):
        with open(f'{a.outdir}/{fn}', 'w') as o:
            o.write(head + '\n')
            for r in rows:
                o.write('\t'.join(r) + '\n')
    write_list('possible_haplotigs.tsv',
               [(g, at.get('possible_haplotig_of', ''), at.get('possible_haplotig', ''), at.get('evidence', ''))
                for g, at in sorted(final_gene_attr.items())
                if 'possible_haplotig_of' in at or 'possible_haplotig' in at],
               'gene\tnear_identical_placed_copy\tstageA_flagged\tevidence')
    for tag, fn in (('possible_gap_fusion', 'gap_fusion_candidates.tsv'),
                    ('possible_gap_misorder', 'gap_misorder_candidates.tsv')):
        write_list(fn, [(fin[t]['gene'], t) for t, nt in sorted(final_tx_note.items()) if tag in nt],
                   'gene\ttranscript')
    rows = []
    for line in open(a.splits):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'gene':
            continue
        rows.append(f + [gene_map.get(f[2], ''), tx_map.get(f[3], '')])
    write_list('splits.tsv', rows, 'refined_gene\told_transcript\tnew_refined_gene\tnew_refined_transcript\t'
               'status\tcut_gaps\tfinal_gene\tfinal_transcript')
    n_pc = sum(1 for g, t in gene_type.items() if t == 'mRNA')
    n_hc = sum(1 for g in hc_genes if gene_type.get(g) == 'mRNA')
    with open(f'{a.outdir}/release_sets_report.tsv', 'w') as o:
        o.write(f'metric\tvalue\nrefined_mrnas\t{len(ref)}\nmapped_mrnas\t{len(tx_map)}\nunmapped_mrnas\t{len(problems)}\n'
                f'protein_coding_genes\t{n_pc}\nhc_protein_coding_genes\t{n_hc}\nhc_proteins\t{n}\n'
                f'mrna_spans_set_to_exon_extent\t{n_span_fixed}\n')
        for t in problems[:50]:
            o.write(f'UNMAPPED\t{t}\n')
    print(open(f'{a.outdir}/release_sets_report.tsv').read())
    if problems:
        sys.exit(f'ERROR: {len(problems)} refined mRNAs not found in the final GFF3')


if __name__ == '__main__':
    main()
