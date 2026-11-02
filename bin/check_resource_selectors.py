#!/usr/bin/env python3
"""Fail if process_resource.config and the modules have drifted apart.

Nextflow only emits a WARNING for a selector that matches no process
("There's no process matching config selector: X"), and `nextflow lint` does not
check config selectors against process names at all. So a renamed or deleted
process leaves a silently dead resource entry -- and, worse, a process with no
selector and no label silently falls back to Nextflow's 1-cpu default.
"""
import glob
import re
import sys

procs = {}
labelled = set()
for path in sorted(glob.glob("modules/*.nf")):
    text = open(path, encoding="utf-8").read()
    for m in re.finditer(r"^process\s+(\w+)\s*\{(.*?)\n\}", text, re.S | re.M):
        procs[m.group(1)] = path
        if re.search(r"label\s+'res_", m.group(2)):
            labelled.add(m.group(1))

cfg = open("process_resource.config", encoding="utf-8").read()
selectors = {s.strip("'\"") for s in re.findall(r"withName:\s*(\S+)\s*\{", cfg)}
labels = {s.strip("'\"") for s in re.findall(r"withLabel:\s*(\S+)\s*\{", cfg)}

problems = []

for sel in sorted(selectors):
    if sel not in procs and not any(ch in sel for ch in ".*'\""):
        problems.append(f"stale withName selector '{sel}' matches no process")

for lab in sorted({l for p in glob.glob('modules/*.nf')
                   for l in re.findall(r"label\s+'(res_[a-z0-9_]+)'",
                                       open(p, encoding='utf-8').read())}):
    if lab not in labels:
        problems.append(f"processes carry label '{lab}' but no withLabel selector defines it")

for proc in sorted(procs):
    if proc not in selectors and proc not in labelled:
        problems.append(
            f"{proc} ({procs[proc]}) has no withName selector and no res_ label "
            f"-- it will silently get Nextflow's 1-cpu default")

if problems:
    print("Resource configuration drift:", file=sys.stderr)
    for p in problems:
        print(f"  - {p}", file=sys.stderr)
    sys.exit(1)

print(f"resource selectors OK: {len(procs)} processes, "
      f"{len(selectors)} withName, {len(labels)} withLabel")
