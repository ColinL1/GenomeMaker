#!/usr/bin/env python3
"""Step 11: release QC gate, checks G1-G11 (master plan section 9.3; Plan 1 adds G11).

Every check writes one or more rows to annotation_qc.tsv: check, level (FAIL/WARN), status
(PASS/FAIL/WARN/NA), value, detail. Exit status 1 if any FAIL-level check fails.
Inputs are given as --key path; checks whose inputs are missing are reported NA and, for FAIL-level
checks, count as failures (the release is not complete without them).

usage: annotation_qc_gate.py --release-gff R.gff3 --refined-gff lifted.gff3 --refine-report refine_report.tsv
   --protein-identity protein_identity.tsv --release-sets-report release_sets_report.tsv
   --final-proteins F.proteins.fa --refined-proteins refined.proteins.faa --id-map id_map.tsv
   --gt-report gt1.txt [--gt-report gt2.txt] --stats-json F.stats.json --summary-tsv summary.tsv
   --readme README.md --multiqc-tsv pipeline_summary_mqc.tsv --name-map scaffold_name_map.tsv
   --agp scaffolds.agp --fasta scaffolds.fa --repeat-lib repeats_classified.lib --relabel-report relabel_report.tsv
   --fcs-adaptor fcs_adaptor_report.txt --fcs-gx x.fcs_gx_report.txt --vecscreen vecscreen_regions.tsv
   [--vecscreen-accepted accepted.tsv] --val-dir ncbi/
   --locus-tag X --sbt path --isolate X [--n-placed 14 --chrom-prefix NC_]
   [--expected-lineage 'Eukaryota; Metazoa' --expected-division INV]
   --out annotation_qc.tsv
"""
import argparse
import bisect
import collections
import glob
import json
import os
import re
import sys

ROWS = []


def row(check, level, ok, value, detail=''):
    status = 'NA' if ok is None else ('PASS' if ok else ('FAIL' if level == 'FAIL' else 'WARN'))
    ROWS.append((check, level, status, str(value), detail))


def attrs(s):
    return dict(kv.split('=', 1) for kv in s.strip().rstrip(';').split(';') if '=' in kv)


def gff_iter(p):
    for line in open(p):
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) >= 9:
            yield f, attrs(f[8])


def fasta(p):
    d, n, b = {}, None, []
    for l in open(p):
        if l.startswith('>'):
            if n:
                d[n] = ''.join(b)
            n, b = l[1:].split()[0], []
        else:
            b.append(l.strip())
    if n:
        d[n] = ''.join(b)
    return d


def have(p):
    return bool(p) and all(x and os.path.exists(x) for x in (p if isinstance(p, list) else [p]))


REVCOMP = str.maketrans('ACGTN', 'TGCAN')
CANONICAL_INTRONS = {('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')}


def model_checks(gff, seqs):
    """G1 checks of the protein-coding models against the released sequence. Returns a dict of
    check -> (number of bad items, examples)."""
    tx, cds, utr, exons = {}, collections.defaultdict(list), collections.defaultdict(list), collections.defaultdict(list)
    for f, at in gff_iter(gff):
        if f[2] == 'mRNA':
            tx[at['ID']] = (f[0], f[6], at['Parent'])
        elif f[2] == 'CDS':
            for p in at['Parent'].split(','):
                cds[p].append((int(f[3]), int(f[4]), int(f[7]) if f[7] in '012' else 0))
        elif f[2] in ('five_prime_UTR', 'three_prime_UTR'):
            for p in at['Parent'].split(','):
                utr[p].append((int(f[3]), int(f[4])))
        elif f[2] == 'exon':
            for p in at['Parent'].split(','):
                exons[p].append((int(f[3]), int(f[4])))
    seqs = {k: v.upper() for k, v in seqs.items()}
    gaps = {k: [(m.start() + 1, m.end()) for m in re.finditer(r'N+', v)] for k, v in seqs.items()}
    bins = collections.defaultdict(list)          # same-strand CDS index, 100 kb bins
    for t, (c, s, g) in tx.items():
        for x, y, _ in cds[t]:
            for b in range(x // 100000, y // 100000 + 1):
                bins[(c, s, b)].append((x, y, g))
    bad = {k: [] for k in ('cds_in_gap', 'partial_end', 'utr_over_other_cds', 'utr_junction')}
    for t, (c, s, g) in tx.items():
        if not cds[t] or c not in seqs:
            continue
        sq, gp = seqs[c], gaps[c]
        for x, y, _ in cds[t]:
            i = bisect.bisect_right(gp, (y, float('inf'))) - 1      # last gap starting at or before y
            if i >= 0 and gp[i][1] >= x:
                bad['cds_in_gap'].append(f'{t}:{x}-{y}')
        for x, y in utr[t]:
            if any(a <= y and b >= x and gg != g for bb in range(x // 100000, y // 100000 + 1)
                   for a, b, gg in bins.get((c, s, bb), ())):
                bad['utr_over_other_cds'].append(f'{t}:{x}-{y}')
        # partial ends as funannotate/NCBI read them: allowed at a sequence end, a gap or a splice consensus
        segs = sorted(cds[t], key=lambda z: z[0], reverse=(s == '-'))
        nt = ''.join(sq[a - 1:b] if s == '+' else sq[a - 1:b][::-1].translate(REVCOMP) for a, b, _ in segs)
        ph = segs[0][2]
        p5, p3 = ph != 0 or nt[:3] != 'ATG', (len(nt) - ph) % 3 != 0 or nt[-3:] not in ('TAA', 'TAG', 'TGA')
        lo, hi = min(x for x, _, _ in cds[t]), max(y for _, y, _ in cds[t])
        for side, p in (('5', p5), ('3', p3)):
            if not p:
                continue
            if (side == '5') == (s == '+'):
                ok = lo == 1 or sq[lo - 2] == 'N' or sq[max(0, lo - 3):lo - 1] == ('AG' if s == '+' else 'AC') or \
                     (s == '-' and sq[max(0, lo - 3):lo - 1] == 'GC')
            else:
                ok = hi == len(sq) or sq[hi] == 'N' or sq[hi:hi + 2] in (('GT', 'GC') if s == '+' else ('CT',))
            if not ok:
                bad['partial_end'].append(f"{t}:{side}'")
        # introns outside the CDS (UTR junctions) must have splice consensus ends
        ex = sorted(exons[t])
        cint = {(b + 1, x - 1) for (a, b, _), (x, y, _) in zip(sorted(cds[t]), sorted(cds[t])[1:])}
        for (a, b), (x, y) in zip(ex, ex[1:]):
            i1, i2 = b + 1, x - 1
            if (i1, i2) in cint:
                continue
            d, acc = sq[i1 - 1:i1 + 1], sq[i2 - 2:i2]
            if s == '-':
                d, acc = acc[::-1].translate(REVCOMP), d[::-1].translate(REVCOMP)
            if (d, acc) not in CANONICAL_INTRONS:
                bad['utr_junction'].append(f'{t}:{i1}-{i2}')
    return bad


def main():
    ap = argparse.ArgumentParser()
    for k in ('release-gff', 'refined-gff', 'refine-report', 'protein-identity', 'release-sets-report',
              'final-proteins', 'refined-proteins', 'id-map', 'stats-json', 'summary-tsv', 'readme',
              'multiqc-tsv', 'name-map', 'agp', 'fasta', 'repeat-lib', 'relabel-report', 'fcs-adaptor',
              'fcs-gx', 'fcs-gx-accepted', 'vecscreen', 'vecscreen-accepted', 'val-dir', 'locus-tag', 'sbt',
              'isolate', 'strain', 'out'):
        ap.add_argument('--' + k)
    ap.add_argument('--gt-report', action='append', default=[])
    ap.add_argument('--screen-actions', help='screen_actions.tsv of the screen that ran before scaffolding')
    ap.add_argument('--screen-review', help='screen_review.tsv of that screen (decisions on the unclear hits)')
    ap.add_argument('--screens-not-run', action='store_true',
                    help='the contamination screens were switched off: G9 is reported as WARN, not run')
    ap.add_argument('--n-placed', type=int, default=0,
                    help='chromosome-scale scaffolds (scaffold_1..N); 0 = do not check their names')
    ap.add_argument('--chrom-prefix', default='',
                    help='name prefix the chromosome scaffolds map to in the name map (e.g. NC_); empty = no check')
    ap.add_argument('--expected-lineage', default='Eukaryota',
                    help='G11: every lineage in the .sqn must start with this')
    ap.add_argument('--expected-division', default='',
                    help='G11: taxonomic division every organism must carry (e.g. INV); empty = no check')
    a = ap.parse_args()

    # ---------------- G1 / G2 on the release GFF
    pc_genes = set()
    if have(a.release_gff):
        mrna, exons, cds = {}, collections.defaultdict(list), collections.defaultdict(list)
        utr = collections.defaultdict(lambda: {'5': [], '3': []})
        agat = 0
        types = collections.Counter()
        trna_iso = collections.Counter()
        ncrna = []
        for f, at in gff_iter(a.release_gff):
            types[f[2]] += 1
            if f[1] == 'AGAT':
                agat += 1
            if f[2] == 'mRNA':
                mrna[at['ID']] = (f[0], f[6], int(f[3]), int(f[4]), at['Parent'])
                pc_genes.add(at['Parent'])
            elif f[2] == 'exon':
                for p in at['Parent'].split(','):
                    exons[p].append((int(f[3]), int(f[4])))
            elif f[2] == 'CDS':
                for p in at['Parent'].split(','):
                    cds[p].append((int(f[3]), int(f[4])))
            elif f[2] in ('five_prime_UTR', 'three_prime_UTR'):
                utr[at['Parent']]['5' if f[2][0] == 'f' else '3'].append(int(f[4]) - int(f[3]) + 1)
            elif f[2] == 'tRNA':
                m = re.search(r'tRNA-([A-Za-z]+)', at.get('product', at.get('Name', '')))
                trna_iso[m.group(1) if m else '?'] += 1
            elif f[2] in ('ncRNA', 'rRNA'):
                ncrna.append((f[0], f[3], f[4], f[6], at))
        row('G1_no_AGAT_features', 'FAIL', agat == 0, agat)
        # partial ends from the refined GFF (start/stop codon features), mapped by transcript
        idmap = {}
        if have(a.id_map):
            for l in open(a.id_map):
                x = l.rstrip('\n').split('\t')
                if x[0] != 'refined_gene' and x[3]:
                    idmap[x[1]] = x[3]
        has_start, has_stop = {}, {}
        if have(a.refined_gff):
            st, sp = set(), set()
            ids = set()
            for f, at in gff_iter(a.refined_gff):
                if f[2] == 'mRNA':
                    ids.add(at['ID'])
                elif f[2] == 'start_codon':
                    st.add(at['Parent'])
                elif f[2] == 'stop_codon':
                    sp.add(at['Parent'])
            for t in ids:
                ft = idmap.get(t)
                if ft:
                    has_start[ft], has_stop[ft] = t in st, t in sp
        bad5 = sum(1 for t, u in utr.items() if u['5'] and t in has_start and not has_start[t])
        bad3 = sum(1 for t, u in utr.items() if u['3'] and t in has_stop and not has_stop[t])
        unmapped_utr = sum(1 for t, u in utr.items() if (u['5'] or u['3']) and t not in has_start)
        row('G1_no_UTR_on_partial_end', 'FAIL', (bad5 + bad3 == 0) and unmapped_utr == 0 and bool(has_start),
            bad5 + bad3, f'5p on no-start {bad5}; 3p on no-stop {bad3}; UTR transcripts not mapped {unmapped_utr}')
        l5 = sorted(x for u in utr.values() for x in u['5'])
        l3 = sorted(x for u in utr.values() for x in u['3'])
        short = sum(1 for x in l5 + l3 if x <= 3)
        row('G1_UTR_lengths', 'FAIL', short == 0,
            f"5p tx {sum(1 for u in utr.values() if u['5'])}; 3p tx {sum(1 for u in utr.values() if u['3'])}",
            f'pieces <=3 bp {short}; median 5p {l5[len(l5) // 2] if l5 else 0}; median 3p {l3[len(l3) // 2] if l3 else 0}')
        span_bad = sum(1 for t, (c, s, a0, b0, g) in mrna.items()
                       if exons[t] and (min(x for x, _ in exons[t]), max(y for _, y in exons[t])) != (a0, b0))
        cds_out = sum(1 for t in mrna for x, y in cds[t] if not any(e0 <= x and y <= e1 for e0, e1 in exons[t]))
        row('G2_transcript_span_equals_exons', 'FAIL', span_bad == 0 and cds_out == 0, span_bad, f'CDS outside exon {cds_out}')
        n_t = sum(trna_iso.values())
        top = trna_iso.most_common(1)[0] if trna_iso else ('-', 0)
        row('G4_tRNA_isotype_share', 'WARN', not n_t or top[1] / n_t <= 0.15, f'{n_t} tRNAs',
            f'max {top[0]} {top[1]} ({100 * top[1] / max(1, n_t):.1f}%)')
        und = sum(v for k, v in trna_iso.items() if k in ('Undet', 'Sup', 'Xxx', '?'))
        row('G4_no_Undet_Sup_tRNA', 'FAIL', und == 0, und)
        row('G4_no_pseudogenes_in_main_GFF', 'FAIL', types['pseudogene'] == 0 and types['RNA'] == 0,
            types['pseudogene'] + types['RNA'])
        mir = sum(1 for x in ncrna if 'miRNA' in x[4].get('ncRNA_class', '') or 'mir-' in x[4].get('product', '').lower()
                  or 'microRNA' in x[4].get('product', ''))
        row('G4_no_miRNA', 'WARN', mir == 0, mir,
            'miRNA calls from Rfam covariance models without small-RNA support; review before submission')
        bycoord = collections.defaultdict(set)
        for c, s, e, strand, _ in ncrna:
            bycoord[(c, s, e)].add(strand)
        dup = sum(1 for v in bycoord.values() if len(v) > 1)
        row('G4_no_opposite_strand_ncRNA', 'FAIL', dup == 0, dup)
        noneuk = sum(1 for x in ncrna if re.search(r'archaea|bacteria|microsporidia', json.dumps(x[4]), re.I))
        row('G4_no_noneukaryotic_rRNA', 'FAIL', noneuk == 0, noneuk)
        row('genes_protein_coding_release_gff', 'INFO', None if not pc_genes else True, len(pc_genes))
    else:
        row('G1_G2_G4_release_gff', 'FAIL', None, 'missing')

    # ---------------- G1 model checks against the released sequence
    if have([a.release_gff, a.fasta]):
        bad = model_checks(a.release_gff, fasta(a.fasta))
        for k, name in (('cds_in_gap', 'G1_no_CDS_in_gap'), ('partial_end', 'G1_partial_ends_valid'),
                        ('utr_over_other_cds', 'G1_no_UTR_over_other_gene_CDS'),
                        ('utr_junction', 'G1_UTR_junctions_canonical')):
            row(name, 'FAIL', not bad[k], len(bad[k]), '; '.join(bad[k][:5]))
    else:
        row('G1_model_checks', 'FAIL', None, 'missing')

    # ---------------- G2 protein identity, G3 refine assertions
    if have(a.protein_identity):
        d = dict(l.rstrip('\n').split('\t')[:2] for l in open(a.protein_identity) if not l.startswith('PROBLEM'))
        row('G2_proteins_identical_except_splits', 'FAIL', d.get('problems') == '0',
            f"identical {d.get('identical_unsplit')}; new {d.get('split_new_proteins')}", f"problems {d.get('problems')}")
    else:
        row('G2_proteins_identical_except_splits', 'FAIL', None, 'missing')
    if have([a.final_proteins, a.refined_proteins, a.id_map]):
        fp, rp = fasta(a.final_proteins), fasta(a.refined_proteins)
        mism = miss = xonly = 0
        for l in open(a.id_map):
            x = l.rstrip('\n').split('\t')
            if x[0] == 'refined_gene':
                continue
            if not x[3] or x[3] not in fp:
                miss += 1
            else:
                # same length, same residues; X (funannotate, codon running into a gap) matches the
                # residue gffread derives from the unambiguous bases of that codon
                a_, b_ = fp[x[3]].rstrip('*.').upper(), rp.get(x[1], '').rstrip('*.').upper()
                if len(a_) != len(b_) or any(p != q and 'X' not in (p, q) for p, q in zip(a_, b_)):
                    mism += 1
                elif a_ != b_:
                    xonly += 1
        row('G2_funannotate_proteins_equal_refined', 'FAIL', mism == 0 and miss == 0, mism,
            f'unmapped {miss}; differing only at X (gap codons) {xonly}')
    if have(a.refine_report):
        fails = [l for l in open(a.refine_report) if l.startswith('FAIL')]
        cnt = {l.split('\t')[1]: l.rstrip('\n').split('\t')[2] for l in open(a.refine_report) if l.startswith('COUNT')}
        row('G3_refine_assertions', 'FAIL', not fails, len(fails))
        row('G3_gene_counts', 'INFO', True, f"BRAKER {cnt.get('braker_genes')} -> AGAT {cnt.get('agat_genes')} -> "
            f"WP2 {cnt.get('genes_after_wp2')} -> WP6 {cnt.get('genes_protein_coding')}",
            f"merges kept {cnt.get('merge_pairs_kept')}, undone {cnt.get('merge_undone')}, added {cnt.get('merge_added')}")
    if a.gt_report:
        bad = [p for p in a.gt_report if not have(p) or re.search(r'warning|error', open(p).read(), re.I)]
        row('GFF3_gt_validator', 'FAIL', not bad, f'{len(a.gt_report) - len(bad)}/{len(a.gt_report)} clean')

    # ---------------- G6 one protein-coding number
    nums = {}
    if have(a.stats_json):
        js = json.load(open(a.stats_json))
        an = js.get('annotation', js)
        g = int(an.get('genes', 0))
        nums['stats.json'] = g - int(an.get('tRNA', 0)) - int(an.get('rRNA', 0)) - int(an.get('ncRNA', 0))
    if pc_genes:
        nums['release_gff'] = len(pc_genes)
    for key, p in (('summary.tsv', a.summary_tsv), ('README', a.readme), ('MultiQC', a.multiqc_tsv)):
        if have(p):
            m = re.search(r'[Pp]rotein-coding genes\D{0,40}?([\d,]+)', open(p).read())
            nums[key] = int(m.group(1).replace(',', '')) if m else None
    vals = set(nums.values())
    row('G6_one_protein_coding_number', 'FAIL', len(nums) == 5 and len(vals) == 1 and None not in vals,
        ';'.join(f'{k}={v}' for k, v in nums.items()))

    # ---------------- G7 names and AGP
    if have(a.name_map):
        st = collections.Counter()
        chrom_ok = 0
        for l in open(a.name_map):
            x = l.rstrip('\n').split('\t')
            if x[0] in ('final_name', '-'):
                continue
            st[x[3]] += 1
            m = re.match(r'scaffold_(\d+)$', x[0])
            if m and int(m.group(1)) <= a.n_placed and x[1].startswith(a.chrom_prefix) and x[3] != 'not_found':
                chrom_ok += 1
        chrom_check = not (a.n_placed and a.chrom_prefix) or chrom_ok == a.n_placed
        row('G7_scaffold_name_map', 'FAIL', st['not_found'] == 0 and chrom_check,
            dict(st), f'chromosome scaffolds mapped to {a.chrom_prefix or "-"}: {chrom_ok}/{a.n_placed}')
    else:
        row('G7_scaffold_name_map', 'FAIL', None, 'missing')
    if have([a.agp, a.fasta]):
        lens = {k: len(v) for k, v in fasta(a.fasta).items()}
        end = collections.Counter()
        gaps_bad = gaps = 0
        for l in open(a.agp):
            if l.startswith('#'):
                continue
            x = l.rstrip('\n').split('\t')
            end[x[0]] = max(end[x[0]], int(x[2]))
            if x[4] in ('N', 'U'):
                gaps += 1
                if not (x[4] == 'U' and x[5] == '100' and x[7] == 'yes' and x[8] == 'align_genus'):
                    gaps_bad += 1
        cover = sum(1 for k, L in lens.items() if end[k] != L)
        row('G7_AGP_gaps_U100_align_genus', 'FAIL', gaps_bad == 0 and cover == 0, f'{gaps} gaps',
            f'bad gap lines {gaps_bad}; scaffolds not covered exactly {cover}')

    # ---------------- G8 repeat library
    if have(a.repeat_lib):
        labs = collections.Counter(l.rstrip('\n').split('#', 1)[1].split()[0] for l in open(a.repeat_lib)
                                   if l.startswith('>') and '#' in l)
        low = sum(v for k, v in labs.items() if re.search(r'(^|/)unknown$', k))
        nonstd = sum(v for k, v in labs.items() if k.startswith(('DNAnona', 'DNAauto')))
        unk = sum(v for k, v in labs.items() if k.split('/')[0] == 'Unknown')
        row('G8_repeat_labels', 'WARN', low == 0 and nonstd == 0, f'lowercase unknown {low}; non-standard {nonstd}',
            f'Unknown families {unk}/{sum(labs.values())}')
    if have(a.relabel_report):
        d = {l.split('\t')[0]: l.rstrip('\n').split('\t')[1:] for l in open(a.relabel_report)}
        row('repeats_unknown_share', 'INFO', True,
            f"before {d.get('unknown_share_of_repeat_bp_before', ['?'])[0]}; after {d.get('unknown_share_of_repeat_bp_after', ['?'])[0]}",
            'share of RepeatMasker-annotated bp labelled Unknown')

    # ---------------- G9 screens on the released FASTA
    def table_rows(p):
        return [l for l in open(p) if l.strip() and not l.startswith('#')]
    screened_once = bool(a.screen_actions or a.screen_review)   # pipeline: one screen, before scaffolding
    if a.screens_not_run:
        for c in ('G9_FCS_adaptor_clean', 'G9_FCS_GX_no_actions', 'G9_VecScreen_no_strong_moderate'):
            ROWS.append((c, 'WARN', 'WARN', 'not run', 'contamination screens switched off (run_screens=false)'))
    elif screened_once:
        if have(a.screen_actions):
            st = collections.Counter(l.split('\t')[0] for l in open(a.screen_actions) if not l.startswith('status'))
            row('G9_screen_actions_applied', 'FAIL', True, f"applied {st['applied']}; listed, not applied {st['not_applied']}",
                'FCS-adaptor, VecScreen, FCS-GX on the contigs before scaffolding; the sequence is not changed afterwards')
        else:
            row('G9_screen_actions_applied', 'FAIL', None, 'missing')
        if have(a.screen_review):
            rv = [l.rstrip('\n').split('\t') for l in open(a.screen_review) if not l.startswith('seq\t')]
            und = sum(1 for x in rv if len(x) > 6 and x[6] == 'UNDECIDED')
            row('G9_screen_hits_decided', 'FAIL', und == 0, f'{len(rv)} hits; undecided {und}',
                '; '.join(f'{k} {v}' for k, v in collections.Counter(x[6] for x in rv if len(x) > 6).items()))
        else:
            row('G9_screen_hits_decided', 'FAIL', None, 'missing')
    elif have(a.fcs_adaptor):
        r = table_rows(a.fcs_adaptor)
        row('G9_FCS_adaptor_clean', 'FAIL', not r, len(r))
    else:
        row('G9_FCS_adaptor_clean', 'FAIL', None, 'missing')
    if a.screens_not_run or screened_once:
        pass
    elif have(a.fcs_gx):
        acc = set()
        if have(a.fcs_gx_accepted):
            acc = {tuple(l.split('\t')[:3]) for l in open(a.fcs_gx_accepted) if not l.startswith('seq')}
        rows_gx = [l.split('\t') for l in table_rows(a.fcs_gx)]
        r = [x for x in rows_gx if x[4] != 'INFO' and tuple(x[:3]) not in acc]
        row('G9_FCS_GX_no_actions', 'FAIL', not r, len(r),
            f'accepted by owner {sum(1 for x in rows_gx if tuple(x[:3]) in acc)}')
    else:
        row('G9_FCS_GX_no_actions', 'FAIL', None, 'missing')
    if a.screens_not_run or screened_once:
        pass
    elif have(a.vecscreen):
        acc = set()
        if have(a.vecscreen_accepted):
            acc = {tuple(l.split('\t')[:3]) for l in open(a.vecscreen_accepted) if not l.startswith('seq')}
        sm = [l.split('\t') for l in open(a.vecscreen) if not l.startswith('seq')]
        bad = [x for x in sm if x[3] in ('Strong', 'Moderate') and tuple(x[:3]) not in acc]
        row('G9_VecScreen_no_strong_moderate', 'FAIL', not bad, len(bad),
            f'accepted by owner {sum(1 for x in sm if tuple(x[:3]) in acc)}; weak {sum(1 for x in sm if x[3] == "Weak")}')
    else:
        row('G9_VecScreen_no_strong_moderate', 'FAIL', None, 'missing')

    # ---------------- G5, G10, G11 NCBI
    g5 = bool(a.locus_tag and have(a.sbt) and (a.isolate or a.strain))
    row('G5_ncbi_metadata', 'FAIL', g5 or None,
        f'locus_tag={a.locus_tag}; sbt={bool(have(a.sbt))}; isolate={a.isolate}; strain={a.strain}')
    if a.val_dir and os.path.isdir(a.val_dir):
        sqn = glob.glob(f'{a.val_dir}/*.sqn')
        lin_bad = div_bad = n_src = 0
        lins = set()
        for p in sqn:
            txt = open(p, errors='replace').read()
            for m in re.finditer(r'orgname\s*\{(.*?)\n\s*\}\s*\}', txt, re.S):
                pass
            for org in re.finditer(r'org\s*\{', txt):
                n_src += 1
            for m in re.finditer(r'lineage\s+"([^"]*)"', txt):
                lins.add(' '.join(m.group(1).split()))
                if not ' '.join(m.group(1).split()).startswith(a.expected_lineage):
                    lin_bad += 1
            # taxonomic division of the organism (orgname: ... gcode N, [mgcode N,] div "X"); other 'div'
            # fields, e.g. the submitter affiliation's department in the .sbt, are not divisions
            divs = re.findall(r'gcode\s+\d+,\s*(?:mgcode\s+\d+,\s*)?(?:pgcode\s+\d+,\s*)?div\s+"([^"]*)"', txt)
            if a.expected_division:
                div_bad += sum(1 for d in divs if d != a.expected_division) + (0 if divs else 1)
            if not re.search(r'lineage\s+"', txt):
                lin_bad += 1
        row('G11_organism_lineage', 'FAIL', bool(sqn) and lin_bad == 0 and div_bad == 0,
            f'{len(sqn)} sqn; lineage bad {lin_bad}; div bad {div_bad}', '; '.join(sorted(lins))[:300])
        val_err = collections.Counter()
        for p in glob.glob(f'{a.val_dir}/*.val'):
            for l in open(p, errors='replace'):
                m = re.match(r'(?i)(error|reject|critical|fatal)\S*: valid \[([A-Za-z0-9_.]+)\]', l)
                if m:
                    val_err[m.group(2)] += 1
        row('G10_val_no_ERROR_REJECT', 'FAIL', bool(glob.glob(f'{a.val_dir}/*.val')) and not val_err,
            sum(val_err.values()), '; '.join(f'{k} {v}' for k, v in val_err.most_common(8)))
        fatal = set()
        for p in glob.glob(f'{a.val_dir}/*.dr'):
            fatal |= {l.strip() for l in open(p, errors='replace') if l.startswith('FATAL')}
        row('G10_dr_no_FATAL', 'FAIL', bool(glob.glob(f'{a.val_dir}/*.dr')) and not fatal, len(fatal), ' | '.join(sorted(fatal)[:5]))
        # warnings and non-fatal discrepancies do not block an NCBI submission: counted for the report only
        val_warn, dr_cat = collections.Counter(), collections.Counter()
        for p in glob.glob(f'{a.val_dir}/*.val'):
            for l in open(p, errors='replace'):
                m = re.match(r'(?i)warning\S*: valid \[([A-Za-z0-9_.]+)\]', l)
                if m:
                    val_warn[m.group(1)] += 1
        DR_INFO = ('COUNT_NUCLEOTIDES', 'FEATURE_COUNT', 'SOURCE_QUALS', 'FEATURE_LIST')
        for p in glob.glob(f'{a.val_dir}/*.dr'):
            for l in open(p, errors='replace'):
                if l.startswith('Detailed Report'):
                    break
                m = re.match(r'([A-Z_]+): (?:(\d+) )?(.*)', l.rstrip('\n'))
                if m and m.group(1) not in DR_INFO and m.group(1) != 'FATAL':
                    dr_cat[f'dr:{m.group(1)}: {m.group(3)}'] += int(m.group(2) or 1)
        # EC numbers table2asn had to split or delete as obsolete (.ecn): the .sqn would differ from the GFF3
        ecn = [l.rstrip('\n') for p in glob.glob(f'{a.val_dir}/*.ecn') for l in open(p, errors='replace') if l.strip()]
        row('G10_ecn_no_EC_changes', 'FAIL', not ecn, len(ecn), '; '.join(ecn[:5]))
        # product names table2asn had to rewrite (.fixedproducts): the .sqn would differ from the GFF3
        fxp = [l.rstrip('\n') for p in glob.glob(f'{a.val_dir}/*.fixedproducts') for l in open(p, errors='replace') if l.strip()]
        row('G10_no_fixedproducts', 'FAIL', not fxp, len(fxp), '; '.join(fxp[:5]))
        row('G10_val_warnings', 'INFO', True, sum(val_warn.values()), '; '.join(f'{k} {v}' for k, v in val_warn.most_common()))
        row('G10_dr_categories', 'INFO', True, sum(dr_cat.values()), '; '.join(f'{k} {v}' for k, v in dr_cat.most_common()))
    else:
        row('G10_G11_ncbi_files', 'FAIL', None, 'missing')

    with open(a.out, 'w') as o:
        o.write('check\tlevel\tstatus\tvalue\tdetail\n')
        for r in ROWS:
            o.write('\t'.join(r) + '\n')
    for r in ROWS:
        print('\t'.join(r))
    failed = [r for r in ROWS if r[1] == 'FAIL' and r[2] in ('FAIL', 'NA')]
    if failed:
        print(f'QC GATE FAILED: {len(failed)} checks', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
