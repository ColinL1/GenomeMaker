#!/usr/bin/env python3
"""Rank assemblies and pick the winner.

Inputs are discovered in the CWD (Nextflow stages them there):
  <name>.gfastats.txt                    -- gfastats report per assembler
  <name>.<lineage>.short_summary.txt     -- BUSCO summary per assembler
  <name>_quast_report/report.tsv         -- QUAST report per assembler (optional)

Environment:
  EXPECTED_GENOME_SIZE  expected haploid size in bp (required: params.genome_size or GenomeScope2)
  SIZE_TOLERANCE        fractional deviation tolerated, e.g. 0.30 for +/-30%
  BUSCO_TIE_BAND        BUSCO completeness difference (pp) treated as a tie, e.g. 2.0

Ranking:
  1. assemblies outside +/-SIZE_TOLERANCE of the expected size are SIZE_OUTLIERs, ranked last;
  2. among the rest, every assembly within BUSCO_TIE_BAND of the best BUSCO C (S+D; duplication
     is left to purge_dups) is a candidate, and the one closest to the expected size wins;
  3. N50 breaks a tie on size.
assembly_ranking.tsv records the size deviation and which rule decided.
"""
import glob
import os
import re
import sys


def parse_busco(path):
    m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\]", open(path).read())
    if not m:
        return None
    return float(m.group(1)), float(m.group(3))  # complete %, duplicated %


def parse_gfastats(path):
    """Return (n50, total_bases). gfastats key names vary slightly by version."""
    text = open(path).read()
    n50 = 0
    for pattern in (r"Scaffold N50\s*:?\s*(\d+)",
                    r"Contig N50\s*:?\s*(\d+)",
                    r"N50\s*:?\s*(\d+)"):
        m = re.search(pattern, text, re.IGNORECASE)
        if m:
            n50 = int(m.group(1))
            break
    total = 0
    for pattern in (r"Total scaffold length\s*:?\s*(\d+)",
                    r"Total contig length\s*:?\s*(\d+)",
                    r"Total length\s*:?\s*(\d+)"):
        m = re.search(pattern, text, re.IGNORECASE)
        if m:
            total = int(m.group(1))
            break
    return n50, total


def parse_quast(name):
    """Pull a few QUAST metrics if its report is present."""
    hits = glob.glob(f"{name}_quast_report/report.tsv") or glob.glob(f"{name}*report.tsv")
    if not hits:
        return {}
    out = {}
    for line in open(hits[0]):
        f = line.rstrip("\n").split("\t")
        if len(f) < 2:
            continue
        key, val = f[0].strip(), f[1].strip()
        if key in ("# contigs", "Largest contig", "GC (%)", "N's per 100 kbp",
                   "Total length", "N50", "L50"):
            try:
                out[key] = float(val)
            except ValueError:
                pass
    return out


expected = float(os.environ.get("EXPECTED_GENOME_SIZE") or 0)
size_tol = float(os.environ.get("SIZE_TOLERANCE") or 0.30)
tie_band = float(os.environ.get("BUSCO_TIE_BAND") or 2.0)
if expected <= 0:
    sys.exit("ERROR: no expected genome size -- set params.genome_size or run GenomeScope2 (run_genomescope)")

busco_by_name = {}
for f in glob.glob("*.short_summary.txt"):
    name = f.replace(".short_summary.txt", "").rsplit(".", 1)[0]
    parsed = parse_busco(f)
    if parsed:
        busco_by_name[name] = parsed

stats_by_name = {}
for f in glob.glob("*.gfastats.txt"):
    stats_by_name[f.replace(".gfastats.txt", "")] = parse_gfastats(f)

if not busco_by_name:
    sys.exit("ERROR: no parseable BUSCO short_summary.txt found -- cannot rank assemblies")

rows = []
for name, (complete_pct, dup_pct) in busco_by_name.items():
    n50, total = stats_by_name.get(name, (0, 0))
    quast = parse_quast(name)
    dev = (total - expected) / expected if total > 0 else None
    if dev is None:
        size_flag = "size_unknown"
    elif abs(dev) > size_tol:
        size_flag = f"SIZE_OUTLIER({dev:+.0%})"
    else:
        size_flag = "ok"
    rows.append(dict(name=name, complete=complete_pct, dup=dup_pct, n50=n50, total=total,
                     dev=dev, size_flag=size_flag, decided_by="",
                     quast_contigs=quast.get("# contigs", ""),
                     quast_largest=quast.get("Largest contig", "")))

# Size outliers are ranked last rather than dropped, so the table still shows them and why.
eligible = [r for r in rows if not r["size_flag"].startswith("SIZE_OUTLIER")] or rows
best_c = max(r["complete"] for r in eligible)
band = [r for r in eligible if best_c - r["complete"] <= tie_band]


def size_key(r):
    return abs(r["dev"]) if r["dev"] is not None else float("inf")


band.sort(key=lambda r: (size_key(r), -r["n50"]))
winner = band[0]
n_out = len(rows) - len(eligible)
rule = f"{n_out} size outlier(s) ranked last; " if n_out else ""
if len(band) == 1:
    winner["decided_by"] = rule + "BUSCO completeness"
elif len(band) > 1 and size_key(band[0]) == size_key(band[1]):
    winner["decided_by"] = rule + f"BUSCO within {tie_band} pp, equal size: N50"
else:
    winner["decided_by"] = rule + f"BUSCO within {tie_band} pp: closest to the expected size"
if len(band) > 1:
    print(f"BUSCO C within {tie_band} pp for {[r['name'] for r in band]}; "
          f"closest to the expected size ({expected:,.0f} bp) -> {winner['name']}")

rest = sorted((r for r in rows if r is not winner),
              key=lambda r: (r["size_flag"].startswith("SIZE_OUTLIER"), r not in band, size_key(r), -r["complete"]))
rows = [winner] + rest

with open("assembly_ranking.tsv", "w") as out:
    out.write("assembler\tcomplete_pct\tduplicated_pct\tN50\ttotal_bp\tsize_dev\t"
              "size_check\tquast_contigs\tquast_largest_contig\tselected\tdecided_by\n")
    for r in rows:
        out.write("\t".join(str(x) for x in (
            r["name"], r["complete"], r["dup"], r["n50"], r["total"],
            f"{r['dev']:+.3f}" if r["dev"] is not None else "NA", r["size_flag"], r["quast_contigs"],
            r["quast_largest"], "yes" if r is winner else "", r["decided_by"])) + "\n")

for r in rows:
    print(f"  {r['name']:<12} C={r['complete']:>6}% D={r['dup']:>5}% "
          f"N50={r['n50']:>12,} size={r['total']:>12,} {r['size_flag']}")

if winner["size_flag"].startswith("SIZE_OUTLIER"):
    sys.stderr.write(
        f"\nWARNING: every candidate deviates from the expected genome size "
        f"({expected:,.0f} bp); winner {winner['name']} is {winner['size_flag']}.\n"
        "Check the GenomeScope estimate and the assemblies before trusting this.\n")

with open("best_assembler.txt", "w") as out:
    out.write(winner["name"])

print(f"\nSelected best assembler: {winner['name']} ({winner['decided_by']})")
