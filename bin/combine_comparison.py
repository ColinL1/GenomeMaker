import glob, json, re, csv

stats = {}

for f in glob.glob("*.gfastats.txt"):
    name = f.split(".gfastats.txt")[0]
    stats.setdefault(name, {})["gfastats"] = {}
    with open(f) as fh:
        for line in fh:
            line = line.strip()
            # gfastats lines look like "Scaffold N50: 12759". Some values
            # (e.g. "Base composition (A:C:G:T): 207962656:...") contain
            # colons of their own, but only the real key/value separator is
            # followed by a space, so split on ": " rather than the first
            # bare colon.
            if ": " not in line:
                continue
            key, _, val = line.partition(": ")
            key, val = key.strip(), val.strip()
            if key and val:
                stats[name]["gfastats"][key] = val

for f in glob.glob("*.short_summary.txt"):
    # filename: <name>.<lineage>.short_summary.txt
    base = f.replace(".short_summary.txt", "")
    name = base.rsplit(".", 1)[0]
    with open(f) as fh:
        text = fh.read()
    m = re.search(
        r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)",
        text,
    )
    stats.setdefault(name, {})["busco"] = {}
    if m:
        stats[name]["busco"] = {
            "complete_pct":   float(m.group(1)),
            "single_pct":     float(m.group(2)),
            "duplicated_pct": float(m.group(3)),
            "fragmented_pct": float(m.group(4)),
            "missing_pct":    float(m.group(5)),
            "n_total":        int(m.group(6)),
        }

all_keys = set()
for d in stats.values():
    all_keys.update(d.get("gfastats", {}).keys())
    all_keys.update(d.get("busco", {}).keys())

names = sorted(stats.keys())

with open("assembly_comparison.tsv", "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["statistic"] + names)
    for key in sorted(all_keys):
        row = [key]
        for name in names:
            val = stats[name].get("gfastats", {}).get(key, stats[name].get("busco", {}).get(key, ""))
            row.append(val)
        writer.writerow(row)

with open("assembly_comparison.json", "w") as f:
    json.dump(stats, f, indent=2)

print(f"Comparison report generated for {len(names)} assemblers")
