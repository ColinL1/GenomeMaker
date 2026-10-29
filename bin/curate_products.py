#!/usr/bin/env python3
"""Product-name curation list from funannotate pass 1.

Genes are selected when
  - their symbol is in the manually checked draft -> product from the draft;
  - their symbol is in Gene2Products.need-curating.txt;
  - their product is a GO term / eggNOG description / truncated name (patterns below), or
  - their product is 'hypothetical protein' while they carry a symbol.
For the last two groups the symbol (funannotate's '_N' suffix removed) is looked up in UniProt
(reviewed entries; human first, else any organism, highest annotation score). The recommended
name becomes the product (action 'fix'); without a hit the gene becomes 'hypothetical protein'
without a symbol (action 'remove'). Lookups are cached in a JSON file.

Genes whose product table2asn had to rewrite (its *.fixedproducts list, optional 6th argument; both
table2asn's 'gnl|ncbi|ID.tN' lines and funannotate tbl2asn's tab lines ending in the gene ID) are
selected as well, so the GFF3/GenBank and the .sqn keep the same product.

Symbols on the need-curating list are looked up in UniProt like the others. The optional draft file
(symbol, product, current product, notes; a missing file means none) gives manually checked products: a gene
with that symbol gets the product, selected or not; with a current product, only genes whose product is
exactly that one.

Products that are exactly a GO term name or an InterPro domain/repeat/site name (go.obo, interpro.tsv of
the funannotate database) are selected as well. When the two files are given they must exist.

With --lookup-all every other gene with a symbol is looked up as well (strategy 'full'): a reviewed
entry gives its recommended name (action 'fix', only when it differs from the current product); without
one the gene keeps funannotate's product and symbol (no row).

Reproducibility: every symbol already in the cache is answered from the cache, so a rerun with the cache
of an earlier run gives the same names. The UniProt release of each lookup (X-UniProt-Release header)
is recorded in the cache under '#uniprot_release'. With --offline no query is made and a symbol
missing from the cache stops the run.

usage: curate_products.py [--offline] [--lookup-all] annotations.txt need-curating.txt product_fix_draft.tsv cache.json out.tsv
       [x.fixedproducts|NONE] [go.obo interpro.tsv]
out.tsv: gene, symbol, old_product, new_product, action (fix|remove), source
"""
import collections
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

OFFLINE = '--offline' in sys.argv
LOOKUP_ALL = '--lookup-all' in sys.argv
argv = [x for x in sys.argv if x not in ('--offline', '--lookup-all')]
ann_p, need_p, draft_p, cache_p, out_p = argv[1:6]
BAD = re.compile(r'(^(positive |negative )?regulation of|activity$|^Belongs to|process$|development$|^response to|'
                 r'binding$|signaling pathway$|^cellular|organization$|assembly$|homeostasis$|differentiation$|'
                 r'transport$|^with |\band$|_|^[^A-Za-z]*$|'
                 # funannotate name-cleaning artefacts: 'gene' replaced by 'protein' inside words,
                 # dropped words, truncation, InterPro domain descriptions, GO phrases
                 r'\w+protein\w+|\bous\b|^Domain (first found|present) in|\.,|\.$|^antigen processing|'
                 r'^[A-Za-z]+(-[A-Za-z]+)+$|identity$|^maintenance of|\([A-Z][a-z]+ [a-z]+\)$|'
                 # funannotate shortenings: 'sub'/'super' for subfamily/superfamily, lost hyphen in 'ATP-dependent'
                 r'\bsub\b|\bsuper\b|(?-i:[A-Z]{2,}dependent)|'
                 # eggNOG sentences and GO-term fragments (NCBI SUSPECT_PRODUCT_NAMES)
                 r'\b(is|are|stimulates|regulates|catalyzes|plays)\b|^ligand for|activity,|acting on|'
                 r'^structural constituent|'
                 # clone / model IDs, stray single letters left by name cleaning, lowercase descriptions
                 r'^[A-Z]{2,5}\d{4,}(-P[A-Z])?$|\bch\d+-\d+|(?-i:\sy\s)|\bue$| -[a-z]|'
                 # InterPro / Pfam domain descriptions instead of protein names
                 r'\b(repeats?|domains?|motifs?)$|, outliers|^(domain|family) (with|of|abundant)|^protein family|'
                 r'^a receptor for|other proteins$|-like proteins$|\bfor \(|^(?-i:[A-Z0-9]{2,8})$|'
                 r'^(?-i:[A-Z]{2,}\d+) protein$|'
                 # eggNOG descriptions led by a lowercase symbol, truncated or lowercase UniProt names
                 r'^(?-i:[a-z]+\d+ )|\bp300 CBP|^regulator of|^uncharacterized protein|substrate in neurons$)', re.I)
# products that are exactly a GO term name or an InterPro domain/repeat/site/superfamily name describe
# a function or a domain, not the protein (optional 7th/8th arguments: go.obo, interpro.tsv from the
# funannotate database; InterPro families are left alone, they are often valid protein names)
TERMS = set()
for p in argv[7:9]:
    if not os.path.exists(p):
        sys.exit(f'ERROR: {p} not found (go.obo / interpro.tsv of the funannotate database)')
if len(argv) > 7:
    TERMS |= {l[6:].strip().lower() for l in open(argv[7]) if l.startswith('name: ')}
if len(argv) > 8:
    for l in open(argv[8]):
        f = l.rstrip('\n').split('\t')
        if len(f) > 2 and f[1] in ('Domain', 'Repeat', 'Homologous_superfamily', 'Conserved_site', 'Active_site',
                                    'Binding_site', 'PTM'):
            TERMS.add(f[2].strip().lower())
fixed_by_table2asn = set()
if len(argv) > 6 and os.path.exists(argv[6]):
    for l in open(argv[6]):
        m = re.search(r'gnl\|ncbi\|([^.\s]+)\.t\d+', l)
        if m:
            fixed_by_table2asn.add(m.group(1))
        elif l.startswith('Changed') and '\t' in l:
            # funannotate's tbl2asn *.fixedproducts: message, old product, location, gene ID
            fixed_by_table2asn.add(l.rstrip('\n').split('\t')[-1].strip())
need = {l.split('\t')[0] for l in (open(need_p) if os.path.exists(need_p) else []) if not l.startswith('#') and l.strip()}
draft = collections.defaultdict(list)    # symbol -> [(product, current product or None)]
for l in (open(draft_p) if os.path.exists(draft_p) else []):
    if not l.startswith('#') and l.strip():
        f = l.rstrip('\n').split('\t') + ['']
        draft[f[0]].append((f[1], f[2] or None))
cache = json.load(open(cache_p)) if os.path.exists(cache_p) else {}


def recommended(names):
    """UniProt recommended name, in NCBI style, as table2asn would rewrite it ('Probable X'/'Uncharacterized X'
    -> 'putative X'; 'ribonucleoprotein protein X' -> 'ribonucleoprotein X')."""
    # UniProt lists alternative names as trailing ' (...)' / ' [...]' groups; parentheses inside the
    # recommended name ('Cap-specific mRNA (nucleoside-2'-O-)-methyltransferase 2') are kept
    n = names.strip()
    while n and n[-1] in ')]':
        depth, i = 0, len(n) - 1
        while i >= 0:                       # opening bracket of the trailing group (any nesting)
            depth += n[i] in ')]'
            depth -= n[i] in '(['
            if depth == 0:
                break
            i -= 1
        if i <= 0 or n[i - 1] != ' ':
            break
        n = n[:i].rstrip()
    n = re.sub(r'\b(\w*protein) protein\b', r'\1', n)
    return re.sub(r'^(Probable|Uncharacterized|Putative) ', 'putative ', n)


def uniprot(sym):
    # entries without 'raw' come from caches written before names were parsed with nested brackets: refetch
    # (offline: used as they are)
    if sym in cache and (cache[sym] is None or 'raw' in cache[sym] or OFFLINE):
        return cache[sym]
    if OFFLINE:
        raise SystemExit(f'ERROR: {sym} is not in the UniProt cache and --offline is set (uniprot_offline)')
    res = None
    for q in (f'gene_exact:"{sym}" AND reviewed:true AND organism_id:9606', f'gene_exact:"{sym}" AND reviewed:true'):
        url = ('https://rest.uniprot.org/uniprotkb/search?format=tsv&size=25&fields=accession,gene_primary,'
               'protein_name,organism_id,annotation_score&query=' + urllib.parse.quote(q))
        for attempt in range(3):
            try:
                with urllib.request.urlopen(url, timeout=60) as r:
                    txt = r.read().decode()
                    rel = r.headers.get('X-UniProt-Release')
                if rel and rel not in cache.setdefault('#uniprot_release', []):
                    cache['#uniprot_release'].append(rel)
                break
            except Exception:
                time.sleep(5)
        else:
            raise SystemExit(f'UniProt lookup failed for {sym}')
        rows = [r.split('\t') for r in txt.strip().split('\n')[1:] if r]
        rows = [r for r in rows if r[1].lower() == sym.lower()]
        if rows:
            rows.sort(key=lambda r: -float(r[4] or 0))
            r = rows[0]
            res = {'acc': r[0], 'name': recommended(r[2]), 'raw': r[2], 'taxid': r[3]}
            break
        time.sleep(0.2)
    cache[sym] = res
    return res


genes = collections.OrderedDict()
hd = None
for l in open(ann_p):
    f = l.rstrip('\n').split('\t')
    if hd is None:
        hd = {k: i for i, k in enumerate(f)}
        continue
    if f[hd['Feature']] != 'mRNA':
        continue
    g = f[hd['GeneID']]
    if g not in genes:
        genes[g] = (f[hd['Name']], re.sub(r', variant \d+$', '', f[hd['Product']]))

cnt = collections.Counter()
with open(out_p, 'w') as o:
    o.write('gene\tsymbol\told_product\tnew_product\taction\tsource\n')
    for g, (n, p) in genes.items():
        base = re.sub(r'_\d+$', '', n)
        own = next((new for new, cur in draft.get(base, []) if cur in (None, p)), None)
        if own is not None:
            o.write(f'{g}\t{n}\t{p}\t{own}\tfix\tdraft (manual check)\n')
            cnt['fix_draft'] += 1
            continue
        if not (base in need or BAD.search(p) or p.lower() in TERMS or (p == 'hypothetical protein' and n)
                or g in fixed_by_table2asn):
            if LOOKUP_ALL and base:
                hit = uniprot(base)
                if not hit:
                    cnt['all_no_entry_kept'] += 1
                elif recommended(hit['name']) == p:
                    cnt['all_same_name'] += 1
                else:
                    o.write(f"{g}\t{n}\t{p}\t{recommended(hit['name'])}\tfix\tUniProt {hit['acc']} (taxid {hit['taxid']}; lookup all)\n")
                    cnt['all_fix_uniprot'] += 1
            continue
        hit = uniprot(base) if base else None
        if hit:
            o.write(f"{g}\t{n}\t{p}\t{recommended(hit['name'])}\tfix\tUniProt {hit['acc']} (taxid {hit['taxid']})\n")
            cnt['fix_uniprot'] += 1
        else:
            o.write(f'{g}\t{n}\t{p}\thypothetical protein\tremove\tno reviewed UniProt entry for symbol\n')
            cnt['remove'] += 1
json.dump(cache, open(cache_p, 'w'), indent=0)
print(dict(cnt))
