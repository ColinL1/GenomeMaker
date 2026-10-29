#!/usr/bin/env python3
"""Build the refined GFF3 (Plan 1 step 3; master plan WP1-WP6, checks G1-G4).

Protein-coding models
  WP1  UTRs from BRAKER's braker_utr.gtf (stringtie2utr) instead of AGAT's: drop UTRs equal to the
       stop codon or <= 3 bp, trim the stop codon off 3'UTRs that start on it, keep a 5'UTR only
       upstream of the CDS and a 3'UTR only downstream of it, and never on a partial end.
       Transcript span = exon extent (removes the AUGUSTUS partial-span extension).
       CDS segments and phases are taken from the validated GFF (what EXTRACT_PROTEINS used) and
       asserted equal to BRAKER CDS + stop codon. 3'-partial CDS that BRAKER4 normalize_gtf cut 3 bp
       short (out-of-frame stop-like end, exon 3 bp longer) get the 3 bp back (cds_3prime_restored).
       5'-partial CDS that AUGUSTUS starts 1 bp inside a scaffold end are taken to the end, phase
       +1 (cds_5prime_to_seq_end). UTR exons beyond the first non-canonical (not GT-AG, GC-AG,
       AT-AC) junction out from the CDS are dropped (utr_dropped_noncanonical_junction).
  WP2  Gene membership recomputed from BRAKER genes: same-strand genes are merged (union-find,
       lowest gene number kept) when two transcripts' CDS overlap with frame agreement >= 0.9.
  WP5  Gene attributes evidence=none|partial|full (gene_support.tsv), possible_haplotig_of=<gene>
       (near-identical copy on a placed scaffold) and possible_haplotig=stageA (>= 50% of the gene
       span in a Stage A HAPLOTIG/OVLP region).
  WP1b UTR exons over the CDS of another gene on the same strand (read-through StringTie transcripts)
       or over an assembly gap are dropped, with every UTR exon beyond them.
  WP6  Junctions with split=yes in the gap tags are cut: every transcript of the gene that spans the
       gap is cut in the gap intron; downstream pieces form gene <gene>_s2 (transcripts <gene>_s2.tN).
       A gene spans its transcripts; its ends are not run out to the gap (a partial end at a splice
       consensus is accepted by NCBI, and a run-out hid the problems below).
       Other tagged junctions give mRNA note=possible_gap_fusion / possible_gap_misorder.
       CDS exons that run through a gap are cut at the gap in the same reading frame. If every
       transcript of the gene does so, the gene is cut there as above; otherwise those transcripts are
       dropped and the ones without a gap in the CDS kept. Cut pieces with a CDS under 90 bp are dropped.
       In each piece, a cut transcript whose CDS lies inside another transcript of the piece (same frame)
       and that has no intron of its own with RNA-seq support (--hints, >= --min-intron-reads) is dropped
       (split_remnant_dropped); a transcript with no exon overlap with the rest of its piece becomes its
       own gene <piece gene>_i<N> (split_island_own_gene).
  Partial ends are read from the sequence (as funannotate does). Transcripts with a partial end inside
       a sequence that is not at a splice consensus are dropped (truncated models), and so are partial
       transcripts with a partial end inside the sequence in genes that have a complete transcript.
       splits.tsv lists every dropped transcript (status dropped).
  Owner edits (--transcript-edits TSV, optional): listed transcripts are dropped, or get a corrected CDS
       (checked: start, whole codons, final stop, no internal stop, canonical introns); applied before
       WP2. splits.tsv status dropped / cds_edited.
ncRNA
  WP3  tRNAs kept only if in the EukHighConfidenceFilter ID list and not Undet/Sup, and given
       product=tRNA-<aa> (funannotate's source for the .tbl product); tRNA pseudogenes are written
       to a separate GFF3.
  WP4  Rfam: drop miRNA (unless --keep-mirna), drop non-eukaryotic SSU/LSU rRNA models, drop
       tRNA/Rfam features overlapping an assembly gap (NCBI FeatureBeginsOrEndsInGap), collapse
       identical-coordinate opposite-strand pairs (keep higher rfam_score), drop ncRNAs > 50% inside
       interspersed repeats (not Simple_repeat/Low_complexity) except rRNA and snRNA.
Output has a ##sequence-region line for every scaffold. Any failed assertion exits non-zero.
"""
import argparse
import collections
import gzip
import re
import sys

REPORT = []          # (category, id, detail)
TRNA_PRODUCT = {'SeC': 'Sec', 'SelCys': 'Sec', 'iMet': 'Met', 'fMet': 'Met'}
COUNTS = collections.Counter()
FAIL = []


def rep(cat, ident, detail=''):
    REPORT.append((cat, ident, detail))
    COUNTS[cat] += 1


def opener(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)


def attrs_of(s):
    return dict(kv.split('=', 1) for kv in s.strip().rstrip(';').split(';') if '=' in kv)


def gnum(g):
    m = re.match(r'g(\d+)', g)
    return int(m.group(1)) if m else 10 ** 9


def merge_iv(ivs):
    """Merge overlapping or touching closed intervals."""
    out = []
    for a, b in sorted(ivs):
        if out and a <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [tuple(x) for x in out]


def subtract_iv(iv, cuts):
    a, b = iv
    pieces = [(a, b)]
    for c, d in cuts:
        nxt = []
        for x, y in pieces:
            if d < x or c > y:
                nxt.append((x, y))
                continue
            if x < c:
                nxt.append((x, c - 1))
            if d < y:
                nxt.append((d + 1, y))
        pieces = nxt
    return pieces


# ----------------------------------------------------------------------------- inputs
def read_lengths(fasta):
    """Scaffold lengths, order, N-run intervals (1-based closed) and upper-case sequences."""
    lens, order, gaps, seqs, name, buf = {}, [], {}, {}, None, []

    def done():
        if name is not None:
            s = ''.join(buf)
            lens[name] = len(s)
            gaps[name] = [(m.start() + 1, m.end()) for m in re.finditer(r'[Nn]+', s)]
            seqs[name] = s.upper()
    for line in open(fasta):
        if line.startswith('>'):
            done()
            name, buf = line[1:].split()[0], []
            order.append(name)
        else:
            buf.append(line.strip())
    done()
    return lens, order, gaps, seqs


def overlaps_gap(gaps, a, b):
    return any(x <= b and y >= a for x, y in gaps)


CANONICAL_INTRONS = {('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')}
REVCOMP = str.maketrans('ACGTN', 'TGCAN')


def canonical_intron(seq, a, b, strand):
    """True if intron a-b (1-based closed) has GT-AG, GC-AG or AT-AC ends on its strand."""
    donor, acceptor = seq[a - 1:a + 1], seq[b - 2:b]
    if strand == '-':
        donor, acceptor = acceptor[::-1].translate(REVCOMP), donor[::-1].translate(REVCOMP)
    return (donor, acceptor) in CANONICAL_INTRONS


def drop_noncanonical_utr(cds, utrs, seq, strand):
    """UTR intervals kept, walking out from the CDS on each side up to the first intron
    (between UTR exons, or between the CDS and a UTR exon) without canonical splice ends.
    Coding introns are not checked here. Returns (kept intervals, first bad intron or None)."""
    lo, hi = min(x for x, _, _ in cds), max(y for _, y, _ in cds)
    ex = merge_iv([(x, y) for x, y, _ in cds] + utrs)
    i_lo = next(i for i, (x, y) in enumerate(ex) if x <= lo <= y)
    i_hi = next(i for i, (x, y) in enumerate(ex) if x <= hi <= y)
    keep_lo, keep_hi, bad = ex[i_lo][0], ex[i_hi][1], None
    for i in range(i_hi, len(ex) - 1):              # right of the CDS
        if not canonical_intron(seq, ex[i][1] + 1, ex[i + 1][0] - 1, strand):
            bad = (ex[i][1] + 1, ex[i + 1][0] - 1)
            break
        keep_hi = ex[i + 1][1]
    for i in range(i_lo, 0, -1):                     # left of the CDS
        if not canonical_intron(seq, ex[i - 1][1] + 1, ex[i][0] - 1, strand):
            bad = bad or (ex[i - 1][1] + 1, ex[i][0] - 1)
            break
        keep_lo = ex[i - 1][0]
    return [u for u in utrs if keep_lo <= u[0] and u[1] <= keep_hi], bad


def drop_utr_outward(cds, utrs, bad):
    """UTR intervals kept, walking out from the CDS on each side up to the first one for which bad(u)
    is true; that one and all beyond it are dropped. Returns (kept, dropped)."""
    lo, hi = min(x for x, _, _ in cds), max(y for _, y, _ in cds)
    kept, dropped = [], []
    for side in (sorted(u for u in utrs if u[0] > hi), sorted((u for u in utrs if u[1] < lo), reverse=True)):
        stop = False
        for u in side:
            stop = stop or bad(u)
            (dropped if stop else kept).append(u)
    return kept, dropped


CDS_BIN = 100000


def cds_index(models, gene_of):
    """(chrom, strand, bin) -> [(a, b, gene)] of all CDS segments."""
    idx = collections.defaultdict(list)
    for t, m in models.items():
        for x, y, _ in m['cds']:
            for bn in range(x // CDS_BIN, y // CDS_BIN + 1):
                idx[(m['chrom'], m['strand'], bn)].append((x, y, gene_of[t]))
    return idx


def hits_other_cds(idx, chrom, strand, gene, a, b):
    return any(x <= b and y >= a and g != gene for bn in range(a // CDS_BIN, b // CDS_BIN + 1)
               for x, y, g in idx.get((chrom, strand, bn), ()))


def cds_gaps(m, gaps):
    """Assembly gaps that a CDS segment of model m runs through."""
    return {gp for gp in gaps for x, y, _ in m['cds'] if gp[0] <= y and gp[1] >= x}


def split_cds_at_gaps(cds, gaps, strand):
    """CDS segments with the gap bases removed. Every piece gets the phase of the original reading
    frame, so the proteins on both sides of a gap are those of the uncut model. Returns (cds, phase
    mismatches of the original segments)."""
    segs = sorted(cds, key=lambda x: x[0], reverse=(strand == '-'))
    phase0, cum, out, bad = segs[0][2], 0, [], 0
    for a, b, ph in segs:
        if (3 - (cum - phase0) % 3) % 3 != ph:
            bad += 1
        for x, y in subtract_iv((a, b), sorted(gaps)):
            off = cum + (x - a if strand == '+' else b - y)
            out.append((x, y, (3 - (off - phase0) % 3) % 3))
        cum += b - a + 1
    return sorted(out), bad


def cut_pieces(m, cuts):
    """Pieces of transcript model m between the gaps cuts (in transcript order):
    [(k, piece model, changed, lo, hi)], lo/hi the piece bounds (+-1e12 when open)."""
    strand = m['strand']
    ext_lo = min([x for x, _, _ in m['cds']] + [u[0] for u in m['utr5'] + m['utr3']])
    ext_hi = max([y for _, y, _ in m['cds']] + [u[1] for u in m['utr5'] + m['utr3']])
    out = []
    for k in range(len(cuts) + 1):
        if strand == '+':
            lo_ = cuts[k - 1][1] + 1 if k > 0 else -10 ** 12
            hi_ = cuts[k][0] - 1 if k < len(cuts) else 10 ** 12
        else:
            hi_ = cuts[k - 1][0] - 1 if k > 0 else 10 ** 12
            lo_ = cuts[k][1] + 1 if k < len(cuts) else -10 ** 12
        cds_k = [(x, y, ph) for x, y, ph in m['cds'] if x >= lo_ and y <= hi_]
        if not cds_k:
            continue
        nm = dict(m)
        nm['cds'] = cds_k
        nm['utr5'] = [u for u in m['utr5'] if u[0] >= lo_ and u[1] <= hi_]
        nm['utr3'] = [u for u in m['utr3'] if u[0] >= lo_ and u[1] <= hi_]
        nm['start'] = [c for c in m['start'] if c[0] >= lo_ and c[1] <= hi_]
        nm['stop'] = [c for c in m['stop'] if c[0] >= lo_ and c[1] <= hi_]
        nm['has_start'] = bool(nm['start'])
        nm['has_stop'] = bool(nm['stop'])
        nm['note'] = list(m['note'])
        nm['cut_lo'] = m['cut_lo'] or lo_ > ext_lo     # the transcript went on beyond this piece
        nm['cut_hi'] = m['cut_hi'] or hi_ < ext_hi
        out.append((k, nm, len(cds_k) != len(m['cds']), lo_, hi_))
    return out


STOPS = ('TAA', 'TAG', 'TGA')


def seq_partial(m, seq):
    """(5' partial, 3' partial) read from the sequence, as funannotate and NCBI decide them:
    complete 5' end = phase 0 and ATG, complete 3' end = whole codons ending in a stop."""
    segs = sorted(m['cds'], key=lambda x: x[0], reverse=(m['strand'] == '-'))
    nt = ''.join(seq[a - 1:b] if m['strand'] == '+' else seq[a - 1:b][::-1].translate(REVCOMP) for a, b, _ in segs)
    ph = segs[0][2]
    return ph != 0 or nt[:3] != 'ATG', (len(nt) - ph) % 3 != 0 or nt[-3:] not in STOPS


def end_site(m, seq, side):
    """Where the 5' or 3' CDS end lies: 'cut' (a split end), 'end' (sequence end or next to a gap),
    'splice' (at a splice consensus, NCBI accepts a partial end there) or 'none'."""
    left = (side == '5') == (m['strand'] == '+')
    if left:
        if m['cut_lo']:
            return 'cut'
        pos = min(x for x, _, _ in m['cds'])
        if pos == 1 or seq[pos - 2] == 'N':
            return 'end'
        di = seq[max(0, pos - 3):pos - 1]
        ok = di == 'AG' if m['strand'] == '+' else di in ('AC', 'GC')
    else:
        if m['cut_hi']:
            return 'cut'
        pos = max(y for _, y, _ in m['cds'])
        if pos == len(seq) or seq[pos] == 'N':
            return 'end'
        di = seq[pos:pos + 2]
        ok = di in ('GT', 'GC') if m['strand'] == '+' else di == 'CT'
    return 'splice' if ok else 'none'


def tnum(t):
    m = re.search(r'\.t(\d+)$', t)
    return int(m.group(1)) if m else 0


def read_braker_utr(path):
    tx = {}
    for line in opener(path):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        if f[2] == 'gene':
            continue
        if f[2] == 'transcript':
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            t = m.group(1) if m else f[8].strip()
            g = re.search(r'gene_id "([^"]+)"', f[8])
            tx[t] = dict(chrom=f[0], strand=f[6], source=f[1], gene=g.group(1) if g else t.rsplit('.', 1)[0],
                         CDS=[], start_codon=[], stop_codon=[], five_prime_UTR=[], three_prime_UTR=[], exon=[])
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m or f[2] not in ('CDS', 'start_codon', 'stop_codon', 'five_prime_UTR', 'three_prime_UTR', 'exon'):
            continue
        tx[m.group(1)][f[2]].append((int(f[3]), int(f[4])))
    return tx


def read_validated(path):
    """Returns (cds, agat_parent, mrna_src, trna_blocks, pseudo_blocks, rfam_blocks)."""
    cds = collections.defaultdict(list)      # mRNA -> [(a, b, phase)]
    parent = {}                              # mRNA -> AGAT gene
    src = {}
    lines_by_id = collections.defaultdict(list)
    children = collections.defaultdict(list)
    tops = []
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        at = attrs_of(f[8])
        if f[2] == 'mRNA':
            parent[at['ID']] = at['Parent']
            src[at['ID']] = f[1]
        elif f[2] == 'CDS' and f[1] in ('AUGUSTUS', 'GeneMark.hmm3', 'gmst'):
            for p in at['Parent'].split(','):
                cds[p].append((int(f[3]), int(f[4]), int(f[7])))
        if f[1] in ('tRNAscan-SE', 'Infernal') or (f[1] == 'AGAT' and f[2] == 'RNA'):
            if 'Parent' in at:
                children[at['Parent']].append(f)
            else:
                tops.append(f)
            if 'ID' in at:
                lines_by_id[at['ID']].append(f)
    trna, pseudo, rfam = [], [], []
    for top in tops:
        tid = attrs_of(top[8])['ID']
        block = [top]
        stack = [tid]
        while stack:
            cur = stack.pop()
            for ch in children.get(cur, []):
                block.append(ch)
                cid = attrs_of(ch[8]).get('ID')
                if cid:
                    stack.append(cid)
        if top[2] == 'pseudogene':
            pseudo.append(block)
        elif top[1] == 'tRNAscan-SE':
            trna.append(block)
        elif top[1] == 'Infernal':
            rfam.append(block)
    return cds, parent, src, trna, pseudo, rfam


def read_rm_intervals(path):
    iv = collections.defaultdict(list)
    for line in open(path):
        f = line.split()
        if len(f) < 11 or not f[0].isdigit():
            continue
        cls = f[10].split('/')[0]
        if cls in ('Simple_repeat', 'Low_complexity'):
            continue
        iv[f[4]].append((int(f[5]), int(f[6])))
    return {k: merge_iv(v) for k, v in iv.items()}


def covered(ivs, a, b):
    return sum(max(0, min(b, y) - max(a, x) + 1) for x, y in ivs if x <= b and y >= a)


# ----------------------------------------------------------------------------- frame logic
def codon_index(cds, strand):
    """[(a, b, k0)] with k0 = codon-frame index of the transcript-5' base of the segment."""
    segs = sorted(cds, key=lambda x: x[0], reverse=(strand == '-'))
    out, k = [], 0
    phase0 = segs[0][2] if segs else 0
    for a, b, _ in segs:
        out.append((a, b, k - phase0))
        k += b - a + 1
    return out


def frame_agreement(c1, c2, strand):
    """Fraction of overlapping CDS bases read in the same codon position."""
    s1, s2 = codon_index(c1, strand), codon_index(c2, strand)
    tot = same = 0
    for a, b, k1 in s1:
        for c, d, k2 in s2:
            lo, hi = max(a, c), min(b, d)
            if lo > hi:
                continue
            n = hi - lo + 1
            if strand == '+':
                d1, d2 = k1 + (lo - a), k2 + (lo - c)
            else:
                d1, d2 = k1 + (b - hi), k2 + (d - hi)
            tot += n
            same += n if (d1 - d2) % 3 == 0 else 0
    return (same / tot) if tot else None, tot


def frame_of_bases(cds, strand):
    """{base: codon position} over the CDS; equal values in two transcripts mean the same reading frame."""
    out = {}
    for a, b, k0 in codon_index(cds, strand):
        for x in range(a, b + 1):
            out[x] = (k0 + (x - a if strand == '+' else b - x)) % 3
    return out


def read_intron_reads(path):
    """(chrom, start, end) -> RNA-seq reads from the BRAKER hints (src=E intron hints, mult)."""
    reads = collections.Counter()
    if not path or path.split('/')[-1].startswith('NO_'):
        return reads
    for line in opener(path):
        f = line.split('\t')
        if len(f) > 8 and f[2] == 'intron' and 'src=E' in f[8]:
            m = re.search(r'mult=(\d+)', f[8])
            reads[(f[0], int(f[3]), int(f[4]))] += int(m.group(1)) if m else 1
    return reads


def model_introns(m):
    ex = merge_iv([(x, y) for x, y, _ in m['cds']] + m['utr5'] + m['utr3'])
    return [(ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1)]


# ----------------------------------------------------------------------------- owner transcript edits
def cds_phases(segs, strand, ph0):
    """[(a, b)] in genomic order -> [(a, b, phase)], phases following the frame from the 5' end."""
    out, ph = {}, ph0
    for a, b in (segs if strand == '+' else segs[::-1]):
        out[(a, b)] = ph
        ph = (3 - (b - a + 1 - ph) % 3) % 3
    return [(a, b, out[(a, b)]) for a, b in segs]


def apply_transcript_edits(path, models, seqs):
    """Owner decisions per transcript (TSV: transcript, action, cds, reason; '#' lines are comments):
    'drop' removes the transcript (its gene must keep another one); 'cds' replaces its CDS segments
    ('a-b,c-d,...', stop codon included) after checking start codon (if the model had one), whole codons,
    a final stop, no internal stop and canonical introns; UTRs no longer outside the new CDS are dropped.
    Returns (splits.tsv rows, dropped transcript IDs)."""
    rows, dropped = [], set()
    if not path or path.split('/')[-1].startswith('NO_'):
        return rows, dropped
    per_gene = collections.Counter(m['bgene'] for m in models.values())
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        f = line.rstrip('\n').split('\t') + ['', '', '']
        t, action, cds_txt, why = f[0].strip(), f[1].strip(), f[2].strip(), f[3].strip()
        if t not in models:
            FAIL.append(f'transcript edit: {t} not among the models')
            continue
        m = models[t]
        if action == 'drop':
            if per_gene[m['bgene']] < 2:
                FAIL.append(f'transcript edit: dropping {t} leaves gene {m["bgene"]} without transcripts')
                continue
            per_gene[m['bgene']] -= 1
            models.pop(t)
            dropped.add(t)
            rep('tx_dropped_owner_edit', t, why)
            rows.append((m['bgene'], t, '', '', 'dropped', 'tx_dropped_owner_edit'))
        elif action == 'cds':
            s, sq = m['strand'], seqs[m['chrom']]
            segs = sorted(tuple(int(x) for x in p.split('-')) for p in cds_txt.split(','))
            ph0 = 0 if m['has_start'] else (m['cds'][0] if s == '+' else m['cds'][-1])[2]
            nt = ''.join(sq[a_ - 1:b_] for a_, b_ in segs)
            if s == '-':
                nt = nt[::-1].translate(REVCOMP)
            codons = [nt[i:i + 3] for i in range(ph0, len(nt), 3)]
            bad = []
            if m['has_start'] and nt[:3] != 'ATG':
                bad.append('no ATG start')
            if (len(nt) - ph0) % 3:
                bad.append('not whole codons')
            if codons[-1] not in STOPS:
                bad.append('no final stop')
            if any(c in STOPS for c in codons[:-1]):
                bad.append('internal stop')
            for (_, b1), (a2, _) in zip(segs, segs[1:]):
                if not canonical_intron(sq, b1 + 1, a2 - 1, s):
                    bad.append(f'intron {b1 + 1}-{a2 - 1} not canonical')
            if bad:
                FAIL.append(f'transcript edit: {t} new CDS rejected: ' + '; '.join(bad))
                continue
            old = ','.join(f'{x}-{y}' for x, y, _ in m['cds'])
            m['cds'] = cds_phases(segs, s, ph0)
            lo, hi = segs[0][0], segs[-1][1]
            m['stop'] = [(hi - 2, hi)] if s == '+' else [(lo, lo + 2)]
            m['has_stop'] = True
            if m['has_start']:
                m['start'] = [(lo, lo + 2)] if s == '+' else [(hi - 2, hi)]
            for kind in ('utr5', 'utr3'):
                up = kind == 'utr5'
                keep = [u for u in m[kind] if (u[1] < lo if (s == '+') == up else u[0] > hi)]
                for u in m[kind]:
                    if u not in keep:
                        rep('utr_dropped_owner_cds_edit', t, f'{kind}:{u[0]}-{u[1]}')
                m[kind] = keep
            rep('cds_replaced_owner_edit', t, f'{old} -> {cds_txt}; {why}')
            rows.append((m['bgene'], t, m['bgene'], t, 'cds_edited', 'cds_replaced_owner_edit'))
        else:
            FAIL.append(f'transcript edit: {t} unknown action {action!r}')
    return rows, dropped


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    for k in ('--braker-utr-gtf', '--validated-gff', '--fasta', '--trna-hc-ids', '--repeatmasker-out',
              '--gene-support', '--haplotig-tags', '--gap-tags', '--out-gff', '--out-pseudo',
              '--out-report', '--out-splits', '--out-tags'):
        ap.add_argument(k, required=True)
    ap.add_argument('--keep-mirna', action='store_true')
    ap.add_argument('--frame-min', type=float, default=0.9)
    ap.add_argument('--repeat-max', type=float, default=0.5)
    ap.add_argument('--haplotig-pid', type=float, default=98)
    ap.add_argument('--haplotig-qcov', type=float, default=0.8)
    ap.add_argument('--stagea-min', type=float, default=0.5)
    ap.add_argument('--transcript-edits', default=None)
    ap.add_argument('--hints', default=None)
    ap.add_argument('--min-intron-reads', type=int, default=3)
    a = ap.parse_args()

    lens, order, all_gaps, seqs = read_lengths(a.fasta)
    btx = read_braker_utr(a.braker_utr_gtf)
    vcds, agat_parent, vsrc, trna_blocks, pseudo_blocks, rfam_blocks = read_validated(a.validated_gff)

    # ---------------- WP1: models
    if set(btx) != set(vcds):
        FAIL.append(f'mRNA ID sets differ: braker {len(btx)}, validated {len(vcds)}, '
                    f'only braker {len(set(btx) - set(vcds))}, only validated {len(set(vcds) - set(btx))}')
    models = {}
    for t, d in btx.items():
        strand = d['strand']
        cds_v = sorted(vcds.get(t, []))
        cds_b = merge_iv(d['CDS'] + d['stop_codon'])
        if merge_iv([(x, y) for x, y, _ in cds_v]) != cds_b:
            FAIL.append(f'{t}: CDS differs from validated GFF')
            rep('assert_cds_mismatch', t)
        has_start, has_stop = bool(d['start_codon']), bool(d['stop_codon'])
        # BRAKER4 normalize_gtf before the frame check trimmed a stop-like triplet off the 3' end of
        # 3'-partial models even out of frame; the exon still ends 3 bp after the CDS. Put the 3 bp back.
        if not has_stop and cds_v:
            ph5 = (cds_v[0] if strand == '+' else cds_v[-1])[2]
            if (sum(y - x + 1 for x, y, _ in cds_v) - ph5) % 3:
                x, y, ph = cds_v[-1] if strand == '+' else cds_v[0]
                if strand == '+' and any(ea <= y and eb == y + 3 for ea, eb in d['exon']):
                    cds_v[-1] = (x, y + 3, ph)
                    rep('cds_3prime_restored', t, f'{x}-{y} -> {x}-{y + 3}')
                elif strand == '-' and any(ea == x - 3 and x <= eb for ea, eb in d['exon']):
                    cds_v[0] = (x - 3, y, ph)
                    rep('cds_3prime_restored', t, f'{x}-{y} -> {x - 3}-{y}')
        # AUGUSTUS starts a 5'-partial model running in from a scaffold end 1 bp inside it
        # (CDS from base 2 on +; to base L-1 on -). NCBI wants a 5' partial at the sequence end:
        # take the base in and raise the phase by one (same reading frame).
        if not has_start and cds_v:
            L, sq = lens[d['chrom']], seqs[d['chrom']]
            if strand == '+' and cds_v[0][0] == 2 and sq[0] != 'N':
                x, y, ph = cds_v[0]
                cds_v[0] = (1, y, (ph + 1) % 3)
                rep('cds_5prime_to_seq_end', t, f'{x}-{y} -> 1-{y}')
            elif strand == '-' and cds_v[-1][1] == L - 1 and sq[-1] != 'N':
                x, y, ph = cds_v[-1]
                cds_v[-1] = (x, L, (ph + 1) % 3)
                rep('cds_5prime_to_seq_end', t, f'{x}-{y} -> {x}-{L}')
        cmin, cmax = min(x for x, _, _ in cds_v), max(y for _, y, _ in cds_v)
        utrs = {'five_prime_UTR': [], 'three_prime_UTR': []}
        for kind in utrs:
            for u in d[kind]:
                if kind == 'three_prime_UTR' and u in d['stop_codon']:
                    COUNTS['utr3_equal_stop_dropped'] += 1
                    continue
                if u[1] - u[0] + 1 <= 3:
                    rep('utr_dropped_le3bp', t, f'{kind}:{u[0]}-{u[1]}')
                    continue
                pieces = subtract_iv(u, [(x, y) for x, y, _ in cds_v])
                if len(pieces) != 1 or pieces[0] != u:
                    if kind == 'three_prime_UTR' and d['stop_codon'] and \
                            sum(y - x + 1 for x, y in pieces) == u[1] - u[0] + 1 - 3 and len(pieces) == 1:
                        COUNTS['utr3_stop_trimmed'] += 1
                    else:
                        rep('utr_dropped_overlaps_cds', t, f'{kind}:{u[0]}-{u[1]}')
                        continue
                for p in pieces:
                    upstream = (p[1] < cmin) if strand == '+' else (p[0] > cmax)
                    downstream = (p[0] > cmax) if strand == '+' else (p[1] < cmin)
                    ok = upstream if kind == 'five_prime_UTR' else downstream
                    if not ok:
                        rep('utr_dropped_wrong_side', t, f'{kind}:{p[0]}-{p[1]}')
                        continue
                    if (kind == 'five_prime_UTR' and not has_start) or (kind == 'three_prime_UTR' and not has_stop):
                        rep('utr_dropped_partial_end', t, f'{kind}:{p[0]}-{p[1]}')
                        continue
                    if p[1] - p[0] + 1 <= 3:
                        rep('utr_dropped_le3bp', t, f'{kind}:{p[0]}-{p[1]} after trim')
                        continue
                    utrs[kind].append(p)
            utrs[kind] = merge_iv(utrs[kind])
        # stringtie2utr can join a UTR exon to the CDS (or to another UTR exon) by an intron
        # RNA-seq does not support, ending in no splice consensus (NCBI NotSpliceConsensus*):
        # keep UTR exons only up to the first non-canonical junction out from the CDS.
        all_utr = utrs['five_prime_UTR'] + utrs['three_prime_UTR']
        if all_utr:
            kept, bad = drop_noncanonical_utr(cds_v, all_utr, seqs[d['chrom']], strand)
            if bad:
                for kind in utrs:
                    for u in utrs[kind]:
                        if u not in kept:
                            rep('utr_dropped_noncanonical_junction', t, f'{kind}:{u[0]}-{u[1]} intron {bad[0]}-{bad[1]}')
                    utrs[kind] = [u for u in utrs[kind] if u in kept]
        models[t] = dict(chrom=d['chrom'], strand=strand, source=d['source'], bgene=d['gene'],
                         cds=cds_v, utr5=utrs['five_prime_UTR'], utr3=utrs['three_prime_UTR'],
                         start=sorted(d['start_codon']), stop=sorted(d['stop_codon']),
                         has_start=has_start, has_stop=has_stop, note=[], split_of=None, cut_lo=False, cut_hi=False)
    # ---------------- owner transcript edits (--transcript-edits): drop a transcript, or replace its CDS
    edit_rows, edit_dropped = apply_transcript_edits(a.transcript_edits, models, seqs)

    COUNTS['tx_with_utr5'] = sum(1 for m in models.values() if m['utr5'])
    COUNTS['tx_with_utr3'] = sum(1 for m in models.values() if m['utr3'])

    # ---------------- WP2: gene membership
    bgenes = collections.defaultdict(list)
    for t, m in models.items():
        bgenes[m['bgene']].append(t)
    env = {}
    for g, ts in bgenes.items():
        m0 = models[ts[0]]
        env[g] = (m0['chrom'], m0['strand'], min(models[t]['cds'][0][0] for t in ts),
                  max(max(y for _, y, _ in models[t]['cds']) for t in ts))
    uf = {g: g for g in bgenes}

    def find(x):
        while uf[x] != x:
            uf[x] = uf[uf[x]]
            x = uf[x]
        return x

    merges = []
    by_cs = collections.defaultdict(list)
    for g, (c, s, lo, hi) in env.items():
        by_cs[(c, s)].append((lo, hi, g))
    for (c, s), lst in by_cs.items():
        lst.sort()
        active = []
        for lo, hi, g in lst:
            active = [x for x in active if x[1] >= lo]
            for lo2, hi2, g2 in active:
                best = None
                for t1 in bgenes[g]:
                    for t2 in bgenes[g2]:
                        fa, n = frame_agreement(models[t1]['cds'], models[t2]['cds'], s)
                        if fa is not None and (best is None or fa > best[0]):
                            best = (fa, n, t1, t2)
                if best and best[0] >= a.frame_min:
                    merges.append((g2, g, best))
                    r1, r2 = find(g), find(g2)
                    if r1 != r2:
                        keep, drop = sorted((r1, r2), key=gnum)
                        uf[drop] = keep
            active.append((lo, hi, g))
    group = collections.defaultdict(list)
    for g in bgenes:
        group[find(g)].append(g)
    final_gene = {}
    for root, members in group.items():
        gid = min(members, key=gnum)
        for g in members:
            for t in bgenes[g]:
                final_gene[t] = gid
    for g2, g, (fa, n, t1, t2) in merges:
        if fa < a.frame_min:
            FAIL.append(f'merge {g2}+{g} frame agreement {fa:.2f}')
    # cross-check against AGAT, as pairs of BRAKER genes placed in one gene
    def gene_pairs(tx2gene):
        grp = collections.defaultdict(set)
        for t, gid in tx2gene.items():
            grp[gid].add(models[t]['bgene'])
        out = set()
        for v in grp.values():
            v = sorted(v, key=gnum)
            out.update((v[i], v[j]) for i in range(len(v)) for j in range(i + 1, len(v)))
        return out, grp
    agat_pairs, agat_groups = gene_pairs({t: p for t, p in agat_parent.items() if t in models})
    mine_pairs, mine_groups = gene_pairs(final_gene)
    best_fa = {}
    for g2, g, (fa, n, t1, t2) in merges:
        best_fa[tuple(sorted((g, g2), key=gnum))] = fa
    for p in sorted(agat_pairs - mine_pairs, key=lambda x: gnum(x[0])):
        b = [models[t]['cds'] for t in bgenes[p[0]]], [models[t]['cds'] for t in bgenes[p[1]]]
        s_ = models[bgenes[p[0]][0]]['strand']
        fas = [frame_agreement(x, y, s_) for x in b[0] for y in b[1]]
        fas = [f for f in fas if f[0] is not None]
        rep('merge_undone', f'{p[0]},{p[1]}',
            f'AGAT merged; best frame agreement {max(fas)[0]:.2f}' if fas else 'AGAT merged; no CDS overlap')
    for p in sorted(mine_pairs - agat_pairs, key=lambda x: gnum(x[0])):
        rep('merge_added', f'{p[0]},{p[1]}', 'frame agreement >= %.2f (direct or transitive); AGAT kept apart' % a.frame_min)
    COUNTS['merge_pairs_agat'] = len(agat_pairs)
    COUNTS['merge_pairs_kept'] = len(agat_pairs & mine_pairs)
    COUNTS['braker_genes'] = len(bgenes)
    COUNTS['agat_genes'] = len(agat_groups)
    COUNTS['genes_after_wp2'] = len(mine_groups)

    # ---------------- WP1b UTRs over another gene's CDS or over a gap
    # stringtie2utr can give two neighbour genes UTRs from one read-through transcript, each running
    # over the other's CDS (NCBI DuplicateFeat when the spans become equal); a UTR exon can also run
    # through an assembly gap. Walking out from the CDS, such a UTR exon and all beyond it are dropped.
    idx = cds_index(models, final_gene)
    for t, m in models.items():
        if not (m['utr5'] or m['utr3']):
            continue
        gp = all_gaps[m['chrom']]

        def why(u, m=m, t=t, gp=gp):
            if overlaps_gap(gp, *u):
                return 'utr_dropped_overlaps_gap'
            if hits_other_cds(idx, m['chrom'], m['strand'], final_gene[t], *u):
                return 'utr_dropped_overlaps_other_gene_cds'
            return None
        kept, dropped = drop_utr_outward(m['cds'], m['utr5'] + m['utr3'], why)
        for u in dropped:
            rep(why(u) or 'utr_dropped_beyond_dropped_utr', t, f'{u[0]}-{u[1]}')
        m['utr5'] = [u for u in m['utr5'] if u in kept]
        m['utr3'] = [u for u in m['utr3'] if u in kept]

    # ---------------- WP5 inputs
    sup = {}
    for line in open(a.gene_support):
        if line.startswith('#') or line.startswith('gene_id'):
            continue
        f = line.rstrip('\n').split('\t')
        ni, isup, ne, esup = int(f[4]), int(f[9]), int(f[11]), int(f[16])
        if isup == 0 and esup == 0:
            lvl = 0
        elif (ni > 0 and isup == ni) or (ni == 0 and esup > 0):
            lvl = 2
        else:
            lvl = 1
        sup[f[1]] = lvl
    hap = {}
    for line in open(a.haplotig_tags):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'gene':
            continue
        hap[f[0]] = (f[3], float(f[4]) if f[4] else 0, float(f[5]) if f[5] else 0, float(f[7]))

    # ---------------- WP6 tags and splits
    split_gaps = collections.defaultdict(set)   # final gene -> {(intron lo, hi)} of passing junctions
    for line in open(a.gap_tags):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'transcript':
            hd = {k: i for i, k in enumerate(f)}
            continue
        t = f[hd['transcript']]
        if t in edit_dropped:
            continue
        if f[hd['split']] == 'yes':
            split_gaps[final_gene[t]].add((int(f[hd['intron_start']]), int(f[hd['intron_end']]), t))
        elif f[hd['tag']]:
            if f[hd['tag']] not in models[t]['note']:
                models[t]['note'].append(f[hd['tag']])

    gaps_by_chrom = all_gaps

    gene_tx = collections.defaultdict(list)
    for t, gid in final_gene.items():
        gene_tx[gid].append(t)
    splits_rows = list(edit_rows)
    ireads = read_intron_reads(a.hints)

    # AUGUSTUS runs some CDS exons through an assembly gap (an X run in the protein, NCBI
    # FeatureCrossesGap). If every transcript of the gene crosses a gap in its CDS, the CDS is cut at the
    # gap, keeping the reading frame on both sides, and the gene is cut there like a WP6 split; otherwise
    # the transcripts without such a gap stay and the gap-crossing ones are dropped.
    exon_cuts = collections.defaultdict(set)   # gene -> gaps it is cut at
    halves = {}                                # transcript -> gaps, dropped (gene keeps gap-free ones)
    for gid, ts in sorted(gene_tx.items(), key=lambda x: (gnum(x[0]), x[0])):
        cg = {t: cds_gaps(models[t], gaps_by_chrom[models[t]['chrom']]) for t in ts}
        if not any(cg.values()):
            continue
        for t in ts:
            if cg[t]:
                models[t]['cds'], nbad = split_cds_at_gaps(models[t]['cds'], cg[t], models[t]['strand'])
                if nbad:
                    FAIL.append(f'{t}: CDS phases do not follow the reading frame ({nbad} segments)')
        if all(cg.values()):
            exon_cuts[gid] = set().union(*cg.values())
        else:
            halves.update({t: cg[t] for t in ts if cg[t]})
    # a gene that keeps a transcript without a gap in its CDS loses the gap-crossing ones: pieces kept
    # beside a complete transcript make the gene partial in mid-sequence (NCBI PartialProblem,
    # GENE_PARTIAL_CONFLICT), owner decision 2026-10-23
    for t, gps in sorted(halves.items(), key=lambda x: (gnum(x[0]), tnum(x[0]))):
        models.pop(t)
        gid = final_gene.pop(t)
        gene_tx[gid].remove(t)
        cut_txt = ';'.join(f'{x}-{y}' for x, y in sorted(gps))
        rep('tx_dropped_gap_in_cds', t, f'gene {gid} keeps transcripts without a gap in the CDS; gaps {cut_txt}')
        splits_rows.append((gid, t, '', '', 'dropped', 'tx_dropped_gap_in_cds'))
    exon_cut_tx = {t for gid in exon_cuts for t in gene_tx[gid]}

    new_models = {}
    for gid in sorted(set(split_gaps) | set(exon_cuts), key=lambda g: (gnum(g), g)):
        chrom = models[gene_tx[gid][0]]['chrom']
        strand = models[gene_tx[gid][0]]['strand']
        cut_gaps = set(exon_cuts.get(gid, ()))
        for lo, hi, t in split_gaps.get(gid, ()):
            inside = [gp for gp in gaps_by_chrom[chrom] if gp[0] >= lo and gp[1] <= hi]
            cut_gaps.update(inside)
        cuts = sorted(cut_gaps, reverse=(strand == '-'))   # transcript order
        # every transcript must place each cut gap in an intron, or entirely outside its span
        ok = True
        for t in gene_tx[gid]:
            m = models[t]
            ex = merge_iv([(x, y) for x, y, _ in m['cds']] + m['utr5'] + m['utr3'])
            for gp in cuts:
                if any(x <= gp[1] and y >= gp[0] for x, y in ex):
                    ok = False
        if not ok:
            rep('split_skipped_gap_in_exon', gid, ';'.join(f'{x}-{y}' for x, y in cuts))
            continue
        cut_txt = ';'.join(f'{x}-{y}' for x, y in cuts)
        pieces = collections.defaultdict(list)  # piece index -> [(old tx, model, cut by the split)]
        gone = {}      # source transcript -> why a piece of it was dropped
        kept = set()   # source transcripts with at least one piece kept
        for t in sorted(gene_tx[gid]):
            for k, nm, changed, lo_, hi_ in cut_pieces(models[t], cuts):
                changed = changed or t in exon_cut_tx
                if changed and sum(y - x + 1 for x, y, _ in nm['cds']) < 90:
                    rep('split_piece_too_short', t, f'piece {k + 1} of {len(cuts) + 1}')
                    gone.setdefault(t, 'split_piece_too_short')
                    continue
                nm['split_of'] = t
                pieces[k].append((t, nm, changed))

        def place(t, nm, changed, g, nt):
            kept.add(t)
            nm = dict(nm)
            nm['gene_override'] = g
            new_models[nt] = nm
            splits_rows.append((gid, t, g, nt, 'protein_changed' if changed else 'moved_unchanged', cut_txt))
        for k in sorted(pieces):
            ng = gid if k == 0 else f'{gid}_s{k + 1}'
            seen, members = {}, []
            for t, nm, changed in pieces[k]:
                key = (tuple(nm['cds']), tuple(nm['utr5']), tuple(nm['utr3']))
                if key in seen:
                    rep('split_duplicate_piece_dropped', t, f'identical to {seen[key]} in {ng}')
                    gone.setdefault(t, 'split_duplicate_piece')
                    continue
                seen[key] = t
                members.append((t, nm, changed))
            # the cut part of a transcript that went on across the gap, inside another transcript of the
            # piece and without an intron of its own with RNA-seq support, adds nothing but a partial end
            frames = {t: frame_of_bases(nm['cds'], strand) for t, nm, _ in members}
            introns = {t: model_introns(nm) for t, nm, _ in members}
            remnant = set()
            for t, nm, changed in members:
                if not changed:
                    continue
                others = [u for u, _, _ in members if u != t and u not in remnant]
                own = [i for i in introns[t] if all(i not in introns[u] for u in others)
                       and ireads[(chrom, *i)] >= a.min_intron_reads]
                inside = [u for u in others if all(frames[u].get(x) == v for x, v in frames[t].items())]
                if inside and not own:
                    remnant.add(t)
                    rep('split_remnant_dropped', t, f'CDS inside {inside[0]} in {ng}; no own intron with RNA-seq support')
                    gone.setdefault(t, 'split_remnant')
            members = [x for x in members if x[0] not in remnant]
            # transcripts without exon overlap form separate groups: the group with the longest CDS stays
            # the piece's gene, every other group becomes its own gene
            exs = {t: merge_iv([(x, y) for x, y, _ in nm['cds']] + nm['utr5'] + nm['utr3']) for t, nm, _ in members}
            grp = {t: t for t, _, _ in members}

            def root(t):
                while grp[t] != t:
                    t = grp[t]
                return t
            for i, (t, _, _) in enumerate(members):
                for u, _, _ in members[i + 1:]:
                    if any(x <= d and c <= y for x, y in exs[t] for c, d in exs[u]):
                        grp[root(u)] = root(t)
            groups = collections.defaultdict(list)
            for item in members:
                groups[root(item[0])].append(item)
            cdslen = lambda item: sum(y - x + 1 for x, y, _ in item[1]['cds'])
            by_len = sorted(groups.values(), key=lambda g_: -max(cdslen(i) for i in g_))
            n = 0
            for t, nm, changed in by_len[0] if by_len else []:
                if k == 0:
                    nt = t
                else:
                    n += 1
                    nt = f'{ng}.t{n}'
                place(t, nm, changed, ng, nt)
            for j, g_ in enumerate(by_len[1:], 1):
                ig = f'{ng}_i{j}'
                rep('split_island_own_gene', ig, ','.join(t for t, _, _ in g_) + f' no exon overlap with the rest of {ng}')
                for n2, (t, nm, changed) in enumerate(g_, 1):
                    place(t, nm, changed, ig, f'{ig}.t{n2}')
        for t in gene_tx[gid]:
            models.pop(t)
            if t not in kept:
                splits_rows.append((gid, t, '', '', 'dropped', gone.get(t, 'split_piece_too_short')))
        rep('split_gene_at_cds_gap' if gid in exon_cuts else 'split', gid,
            f'{len(pieces)} pieces at gaps ' + ';'.join(f'{x}-{y}' for x, y in cuts))
    for nt, nm in new_models.items():
        models[nt] = nm
        final_gene[nt] = nm['gene_override']
    final_gene = {t: g for t, g in final_gene.items() if t in models}

    # ---------------- partial CDS ends inside a sequence (NCBI PartialProblem*)
    # Partial ends are read from the sequence, as funannotate sets them. A transcript whose partial end
    # is neither at a sequence end, a gap, a cut nor a splice consensus is a truncated model: dropped.
    # In a gene with a complete transcript, partial transcripts whose partial end lies inside the
    # sequence are dropped too (they make the gene partial in the middle of the sequence).
    # protein-coding models with a CDS exon overlapping an rRNA gene are spurious calls in the rDNA (NCBI RNA_CDS_OVERLAP)
    rrna_iv = collections.defaultdict(list)
    for block in rfam_blocks:
        for f in block:
            if f[2] == 'rRNA' and not re.search(r'_(archaea|bacteria|microsporidia)$', attrs_of(f[8]).get('Name', '')):
                rrna_iv[f[0]].append((int(f[3]), int(f[4])))
    for t in sorted(models):
        m = models[t]
        hit = [(x, y) for x, y, _ in m['cds'] if any(a <= y and x <= b for a, b in rrna_iv.get(m['chrom'], ()))]
        if hit:
            rep('tx_dropped_cds_in_rrna', t, 'CDS ' + ','.join(f'{x}-{y}' for x, y in hit))
            g = final_gene.pop(t)
            models.pop(t)
            splits_rows.append((g, t, '', '', 'dropped', 'tx_dropped_cds_in_rrna'))
            if g not in final_gene.values():
                rep('gene_dropped_no_transcripts', g)

    part = {t: seq_partial(m, seqs[m['chrom']]) for t, m in models.items()}
    sites = {t: [end_site(models[t], seqs[models[t]['chrom']], sd) for sd, p in zip('53', part[t]) if p] for t in models}
    gene_tx = collections.defaultdict(list)
    for t, g in final_gene.items():
        gene_tx[g].append(t)
    for g, ts in sorted(gene_tx.items(), key=lambda x: (gnum(x[0]), x[0])):
        complete = [t for t in ts if not any(part[t])]
        for t in sorted(ts):
            if 'none' in sites[t]:
                rep('tx_dropped_truncated', t, 'partial end not at a sequence end, gap or splice consensus')
            elif complete and 'splice' in sites[t]:
                rep('tx_dropped_partial_isoform', t, f'partial end inside the sequence; complete: {",".join(sorted(complete))}')
            else:
                continue
            models.pop(t)
            final_gene.pop(t)
            splits_rows.append((g, t, '', '', 'dropped', REPORT[-1][0]))
        if not any(t in models for t in ts):
            rep('gene_dropped_no_transcripts', g)

    # ---------------- assertions on models (G1, G2)
    for t, m in models.items():
        L = lens[m['chrom']]
        ex = merge_iv([(x, y) for x, y, _ in m['cds']] + m['utr5'] + m['utr3'])
        m['exons'] = ex
        if ex[0][0] < 1 or ex[-1][1] > L:
            FAIL.append(f'{t}: feature beyond scaffold end')
        for x, y, _ in m['cds']:
            if not any(ea <= x and y <= eb for ea, eb in ex):
                FAIL.append(f'{t}: CDS not inside an exon')
        if m['utr5'] and not m['has_start']:
            FAIL.append(f'{t}: 5pUTR on transcript without start codon')
        if m['utr3'] and not m['has_stop']:
            FAIL.append(f'{t}: 3pUTR on transcript without stop codon')
        for u in m['utr3']:
            if u in m['stop']:
                FAIL.append(f'{t}: 3pUTR equals stop codon')

    # checks of the rules above on the final models
    idx = cds_index(models, final_gene)
    for t, m in models.items():
        gp = all_gaps[m['chrom']]
        if any(overlaps_gap(gp, x, y) for x, y, _ in m['cds']):
            FAIL.append(f'{t}: CDS runs through a gap')
        for u in m['utr5'] + m['utr3']:
            if overlaps_gap(gp, *u):
                FAIL.append(f'{t}: UTR runs through a gap')
            if hits_other_cds(idx, m['chrom'], m['strand'], final_gene[t], *u):
                FAIL.append(f'{t}: UTR over the CDS of another gene')
        p = seq_partial(m, seqs[m['chrom']])
        if any(p) and 'none' in [end_site(m, seqs[m['chrom']], sd) for sd, q in zip('53', p) if q]:
            FAIL.append(f'{t}: partial end not at a sequence end, gap or splice consensus')

    # ---------------- gene-level tags
    gene_members = collections.defaultdict(list)
    for t, g in final_gene.items():
        gene_members[g].append(t)
    gene_attr = {}
    bgene_to_final = {}
    for t, g in final_gene.items():
        bgene_to_final.setdefault(models[t]['bgene'], g)
    for g, ts in gene_members.items():
        lv = max(sup.get(models[t]['split_of'] or t, 0) for t in ts)
        at = {'evidence': ['none', 'partial', 'full'][lv]}
        for bg in sorted({models[t]['bgene'] for t in ts}, key=gnum):
            h = hap.get(bg)
            if not h:
                continue
            if h[0] and h[1] >= a.haplotig_pid and h[2] >= a.haplotig_qcov and 'possible_haplotig_of' not in at:
                at['possible_haplotig_of'] = bgene_to_final.get(h[0], h[0])
            if h[3] >= a.stagea_min:
                at['possible_haplotig'] = 'stageA'
        gene_attr[g] = at
    COUNTS['genes_evidence_none'] = sum(1 for v in gene_attr.values() if v['evidence'] == 'none')
    COUNTS['genes_possible_haplotig_of'] = sum(1 for v in gene_attr.values() if 'possible_haplotig_of' in v)
    COUNTS['genes_possible_haplotig_stageA'] = sum(1 for v in gene_attr.values() if 'possible_haplotig' in v)
    COUNTS['mrna_possible_gap_fusion'] = sum(1 for m in models.values() if 'possible_gap_fusion' in m['note'])
    COUNTS['mrna_possible_gap_misorder'] = sum(1 for m in models.values() if 'possible_gap_misorder' in m['note'])

    # ---------------- WP3 tRNAs
    hc = {l.strip() for l in open(a.trna_hc_ids) if l.strip()}
    keep_trna = []
    iso = collections.Counter()
    for block in trna_blocks:
        rna = [f for f in block if f[2] == 'tRNA']
        if not rna:
            FAIL.append('tRNA block without tRNA feature')
            continue
        at = attrs_of(rna[0][8])
        tid, isotype = at['ID'], at.get('isotype', '')
        if tid not in hc:
            rep('trna_dropped_not_hc', tid, isotype)
            continue
        if isotype in ('Undet', 'Sup', 'Supres'):
            rep('trna_dropped_undet_sup', tid, isotype)
            continue
        if overlaps_gap(all_gaps[rna[0][0]], int(rna[0][3]), int(rna[0][4])):
            rep('trna_dropped_overlaps_gap', tid, isotype)
            continue
        iso[isotype] += 1
        # funannotate takes the tRNA product from product=; without it the .tbl says 'product None'
        prod = 'tRNA-' + TRNA_PRODUCT.get(isotype, isotype)
        block = [f[:8] + [f[8].rstrip(';') + f';product={prod}'] if f[2] == 'tRNA' else f for f in block]
        keep_trna.append(block)
    COUNTS['trna_kept'] = len(keep_trna)
    for i, n in iso.items():
        if n / max(1, len(keep_trna)) > 0.15:
            rep('WARN_trna_isotype_gt15pct', i, f'{n}/{len(keep_trna)}')
    COUNTS['trna_pseudogenes_moved'] = len(pseudo_blocks)

    # ---------------- WP4 Rfam
    rm = read_rm_intervals(a.repeatmasker_out)
    rf = []
    for block in rfam_blocks:
        feat = [f for f in block if f[2] in ('ncRNA', 'rRNA')]
        if not feat:
            continue
        at = attrs_of(feat[0][8])
        name, cls = at.get('Name', ''), at.get('ncRNA_class', '')
        if feat[0][2] == 'ncRNA' and cls == 'miRNA' and not a.keep_mirna:
            rep('rfam_dropped_mirna', at['ID'], name)
            continue
        if feat[0][2] == 'rRNA' and re.search(r'_(archaea|bacteria|microsporidia)$', name):
            rep('rfam_dropped_noneuk_rrna', at['ID'], name)
            continue
        if overlaps_gap(all_gaps[feat[0][0]], int(feat[0][3]), int(feat[0][4])):
            rep('rfam_dropped_overlaps_gap', at['ID'], name)
            continue
        rf.append((block, feat[0], at))
    bycoord = collections.defaultdict(list)
    for item in rf:
        f = item[1]
        bycoord[(f[0], f[3], f[4])].append(item)
    rf2 = []
    for key, items in bycoord.items():
        if len(items) > 1 and len({i[1][6] for i in items}) > 1:
            items.sort(key=lambda i: -float(i[2].get('rfam_score', 0)))
            for i in items[1:]:
                rep('rfam_dropped_opposite_strand_dup', i[2]['ID'], f'kept {items[0][2]["ID"]}')
            rf2.append(items[0])
        else:
            rf2.extend(items)
    keep_rfam = []
    for block, f, at in rf2:
        is_rrna = f[2] == 'rRNA'
        is_snrna = at.get('ncRNA_class', '') == 'snRNA'
        a0, b0 = int(f[3]), int(f[4])
        frac = covered(rm.get(f[0], []), a0, b0) / (b0 - a0 + 1)
        if frac > a.repeat_max and not (is_rrna or is_snrna):
            rep('rfam_dropped_in_repeat', at['ID'], f'{at.get("Name", "")} {frac:.2f}')
            continue
        keep_rfam.append(block)
    COUNTS['rfam_kept'] = len(keep_rfam)
    # G4: no opposite-strand duplicates remain
    seen = collections.defaultdict(set)
    for block in keep_rfam:
        f = [x for x in block if x[2] in ('ncRNA', 'rRNA')][0]
        seen[(f[0], f[3], f[4])].add(f[6])
    if any(len(v) > 1 for v in seen.values()):
        FAIL.append('opposite-strand duplicate ncRNAs remain')

    # ---------------- write GFF3
    recs = collections.defaultdict(list)  # chrom -> [(start, sort key, lines)]
    for g, ts in gene_members.items():
        ts = sorted(ts, key=lambda t: (int(re.search(r'\.t(\d+)$', t).group(1)) if re.search(r'\.t(\d+)$', t) else 0))
        m0 = models[ts[0]]
        gs = min(models[t]['exons'][0][0] for t in ts)
        ge = max(models[t]['exons'][-1][1] for t in ts)
        c, s = m0['chrom'], m0['strand']
        ga = ';'.join([f'ID={g}'] + [f'{k}={v}' for k, v in gene_attr[g].items()])
        lines = [f'{c}\t{m0["source"]}\tgene\t{gs}\t{ge}\t.\t{s}\t.\t{ga}']
        for t in ts:
            m = models[t]
            ex = m['exons']
            part = []
            if not m['has_start']:
                part.append('5prime')
            if not m['has_stop']:
                part.append('3prime')
            ma = [f'ID={t}', f'Parent={g}']
            if part:
                ma.append('partial=' + ('both' if len(part) == 2 else part[0]))
            if m['note']:
                ma.append('note=' + ','.join(m['note']))
            lines.append(f'{c}\t{m["source"]}\tmRNA\t{ex[0][0]}\t{ex[-1][1]}\t.\t{s}\t.\t' + ';'.join(ma))
            exo = ex if s == '+' else ex[::-1]
            for k, (x, y) in enumerate(exo, 1):
                lines.append(f'{c}\t{m["source"]}\texon\t{x}\t{y}\t.\t{s}\t.\tID={t}.exon{k};Parent={t}')
            u5 = m['utr5'] if s == '+' else m['utr5'][::-1]
            for k, (x, y) in enumerate(u5, 1):
                lines.append(f'{c}\tstringtie2utr\tfive_prime_UTR\t{x}\t{y}\t.\t{s}\t.\tID={t}.utr5p{k};Parent={t}')
            for x, y in m['start']:
                lines.append(f'{c}\t{m["source"]}\tstart_codon\t{x}\t{y}\t.\t{s}\t0\tParent={t}')
            cd = sorted(m['cds'], reverse=(s == '-'))
            for x, y, ph in cd:
                lines.append(f'{c}\t{m["source"]}\tCDS\t{x}\t{y}\t.\t{s}\t{ph}\tID={t}.cds;Parent={t}')
            for x, y in m['stop']:
                lines.append(f'{c}\t{m["source"]}\tstop_codon\t{x}\t{y}\t.\t{s}\t0\tParent={t}')
            u3 = m['utr3'] if s == '+' else m['utr3'][::-1]
            for k, (x, y) in enumerate(u3, 1):
                lines.append(f'{c}\tstringtie2utr\tthree_prime_UTR\t{x}\t{y}\t.\t{s}\t.\tID={t}.utr3p{k};Parent={t}')
        recs[c].append((gs, ge, g, lines))
    for block in keep_trna + keep_rfam:
        top = block[0]
        recs[top[0]].append((int(top[3]), int(top[4]), attrs_of(top[8])['ID'], ['\t'.join(f) for f in block]))
    for c, lst in recs.items():
        L = lens[c]
        for gs, ge, gid, lines in lst:
            for ln in lines:
                f = ln.split('\t')
                if int(f[3]) < 1 or int(f[4]) > L:
                    FAIL.append(f'{gid}: feature beyond scaffold end')
    with open(a.out_gff, 'w') as o:
        o.write('##gff-version 3\n')
        for c in order:
            o.write(f'##sequence-region {c} 1 {lens[c]}\n')
        for c in order:
            for gs, ge, gid, lines in sorted(recs.get(c, []), key=lambda r: (r[0], -r[1], r[2])):
                o.write('\n'.join(lines) + '\n')
    with open(a.out_pseudo, 'w') as o:
        o.write('##gff-version 3\n')
        used = sorted({b[0][0] for b in pseudo_blocks}, key=order.index)
        for c in used:
            o.write(f'##sequence-region {c} 1 {lens[c]}\n')
        for b in sorted(pseudo_blocks, key=lambda b: (order.index(b[0][0]), int(b[0][3]))):
            o.write('\n'.join('\t'.join(f) for f in b) + '\n')
    with open(a.out_splits, 'w') as o:
        o.write('gene\told_transcript\tnew_gene\tnew_transcript\tstatus\tcut_gaps\n')
        for r in splits_rows:
            o.write('\t'.join(r) + '\n')
    with open(a.out_tags, 'w') as o:
        o.write('gene\ttranscript\tevidence\tpossible_haplotig_of\tpossible_haplotig\tnote\n')
        for g, ts in sorted(gene_members.items(), key=lambda x: (gnum(x[0]), x[0])):
            at = gene_attr[g]
            for t in sorted(ts):
                o.write(f"{g}\t{t}\t{at['evidence']}\t{at.get('possible_haplotig_of', '')}\t"
                        f"{at.get('possible_haplotig', '')}\t{','.join(models[t]['note'])}\n")

    COUNTS['genes_protein_coding'] = len(gene_members)
    COUNTS['mrnas'] = len(models)
    COUNTS['tx_with_utr5_final'] = sum(1 for m in models.values() if m['utr5'])
    COUNTS['tx_with_utr3_final'] = sum(1 for m in models.values() if m['utr3'])
    COUNTS['split_transcript_rows'] = len(splits_rows)
    with open(a.out_report, 'w') as o:
        o.write('category\tid\tdetail\n')
        for k in sorted(COUNTS):
            o.write(f'COUNT\t{k}\t{COUNTS[k]}\n')
        for f in FAIL:
            o.write(f'FAIL\t\t{f}\n')
        for r in REPORT:
            o.write('\t'.join(str(x) for x in r) + '\n')
    for k in sorted(COUNTS):
        print(f'{k}\t{COUNTS[k]}')
    if FAIL:
        print(f'{len(FAIL)} assertion failures, first: {FAIL[:5]}', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
