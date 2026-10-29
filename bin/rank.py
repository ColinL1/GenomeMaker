"""Rank assemblies and pick the winner.

Inputs are discovered in the CWD (Nextflow stages them there):
  <name>.gfastats.txt                    -- gfastats report per assembler
  <name>.<lineage>.short_summary.txt     -- BUSCO summary per assembler
  <name>_quast_report/report.tsv         -- QUAST report per assembler (optional)

Environment:
  EXPECTED_GENOME_SIZE  expected haploid size in bp (0/unset disables the check)
  SIZE_TOLERANCE        fractional deviation tolerated, e.g. 0.30 for +/-30%
  BUSCO_TIE_BAND        BUSCO score difference (pp) treated as a tie, e.g. 0.5

Scoring: BUSCO completeness penalised for duplication, with N50 as a genuine
tie-break inside a tolerance band. Previously N50 only broke EXACT float ties,
which essentially never happen, so contiguity never influenced the choice.
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
    """Pull a few QUAST metrics if its report is present. QUAST was previously
    computed for every assembler and then used by nothing at all."""
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
tie_band = float(os.environ.get("BUSCO_TIE_BAND") or 0.5)

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
    # Primary score rewards completeness and penalises duplication (unpurged
    # haplotigs). N50 is NOT folded in here -- it is applied as a tie-break below,
    # so a marginally-more-complete but far more fragmented assembly cannot win on
    # a rounding difference.
    score = complete_pct - 2 * dup_pct

    size_flag = "ok"
    if expected > 0 and total > 0:
        dev = abs(total - expected) / expected
        if dev > size_tol:
            size_flag = f"SIZE_OUTLIER(+{dev:.0%})" if total > expected else f"SIZE_OUTLIER(-{dev:.0%})"
    elif expected > 0:
        size_flag = "size_unknown"

    rows.append(dict(name=name, complete=complete_pct, dup=dup_pct, n50=n50,
                     total=total, score=score, size_flag=size_flag,
                     quast_contigs=quast.get("# contigs", ""),
                     quast_largest=quast.get("Largest contig", "")))

# Assemblies whose size is wildly off the expectation are ranked last rather than
# silently dropped, so the ranking table still shows them and why.
def sort_key(r):
    return (0 if r["size_flag"] in ("ok", "size_unknown") else 1,
            -r["score"], -r["n50"])

rows.sort(key=sort_key)

# --- tie-break: within tie_band percentage points of the best score, prefer N50 ---
eligible = [r for r in rows if r["size_flag"] in ("ok", "size_unknown")] or rows
best_score = eligible[0]["score"]
tied = [r for r in eligible if best_score - r["score"] <= tie_band]
tied.sort(key=lambda r: -r["n50"])
winner = tied[0]
if len(tied) > 1:
    print(f"BUSCO scores within {tie_band} pp for {[r['name'] for r in tied]}; "
          f"broke the tie on N50 -> {winner['name']}")

with open("assembly_ranking.tsv", "w") as out:
    out.write("assembler\tcomplete_pct\tduplicated_pct\tN50\ttotal_bp\tscore\t"
              "size_check\tquast_contigs\tquast_largest_contig\tselected\n")
    for r in rows:
        out.write("\t".join(str(x) for x in (
            r["name"], r["complete"], r["dup"], r["n50"], r["total"],
            f"{r['score']:.2f}", r["size_flag"], r["quast_contigs"],
            r["quast_largest"], "yes" if r is winner else "")) + "\n")

for r in rows:
    print(f"  {r['name']:<12} C={r['complete']:>6}% D={r['dup']:>5}% "
          f"N50={r['n50']:>12,} size={r['total']:>12,} score={r['score']:>7.2f} {r['size_flag']}")

if winner["size_flag"].startswith("SIZE_OUTLIER"):
    sys.stderr.write(
        f"\nWARNING: every candidate deviates from the expected genome size "
        f"({expected:,.0f} bp); winner {winner['name']} is {winner['size_flag']}.\n"
        "Check the GenomeScope estimate and the assemblies before trusting this.\n")

with open("best_assembler.txt", "w") as out:
    out.write(winner["name"])

print(f"\nSelected best assembler: {winner['name']}")
