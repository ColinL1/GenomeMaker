#!/usr/bin/env python3
"""Methods section of the release README with the tool versions this run recorded.

The template (written by the RELEASE subworkflow from the steps and settings the run used) names
tools as {{tool}}, the key of the tool in the processes' versions.yml files. Each is replaced by the
version(s) the run's processes reported ('a / b' when processes disagree); a tool no process reported
becomes '(version not recorded)'.

usage: release_methods.py methods_template.md all_versions.yml methods.md
"""
import collections
import re
import sys

tpl_p, ver_p, out_p = sys.argv[1:4]

vers = collections.defaultdict(set)
for l in open(ver_p, errors='replace'):
    m = re.match(r'^\s+([^:\s]+):\s*"?([^"]*?)"?\s*$', l.rstrip('\n'))
    if m and m.group(1) != 'container' and m.group(2) and m.group(2).lower() not in ('unknown', 'stub'):
        vers[m.group(1)].add(m.group(2))


def fill(m):
    v = vers.get(m.group(1))
    return ' / '.join(sorted(v)) if v else '(version not recorded)'


txt = re.sub(r'\{\{([^}]+)\}\}', fill, open(tpl_p).read())
open(out_p, 'w').write(txt)
print(txt)
