#!/usr/bin/env python3
"""Prepare a blind whole-course audit: reproducible stratified 20% sample, blind packets, TASK.md.

Usage (from the repository root): python3 private/audits/tools/prepare_course_audit.py COURSE SEED OUTDIR
The population is every ai-draft unit whose course or also_courses is COURSE. Strata are
example-sheet units and Tripos units; each stratum is sorted by id and sampled with
random.Random(SEED).sample(stratum, ceil(0.2 * len(stratum))), sheet stratum first.
"""
import glob, math, os, random, re, sys
import yaml

course, seed, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
units = []
for m in sorted(glob.glob('courses/**/meta.yaml', recursive=True) + glob.glob('tripos/**/meta.yaml', recursive=True)):
    if '/history/' in m or '/variants/' in m:
        continue
    d = yaml.safe_load(open(m))
    if d['solution']['status'] != 'ai-draft':
        continue
    if d['course'] == course or course in (d.get('also_courses') or []):
        units.append(dict(id=d['id'], type=d['type'], dir=os.path.dirname(m), src=d['source']))
if os.path.exists(out):
    sys.exit(f'{out} exists; refusing to redraw an existing sample')
rng = random.Random(seed)
sample, strata = [], {}
for t in ('sheet', 'tripos'):
    pop = sorted([u for u in units if u['type'] == t], key=lambda u: u['id'])
    k = math.ceil(0.2 * len(pop)) if pop else 0
    strata[t] = {'population': len(pop), 'sampled': k}
    sample += rng.sample(pop, k)
os.makedirs(out + '/packets'); os.makedirs(out + '/rework')
rows = []
for u in sample:
    slug = u['id'].replace('/', '_')
    q = open(u['dir'] + '/question.tex').read()
    figs = sorted(set(re.findall(r'\\diagram\{([A-Za-z0-9_-]+)\}', q)))
    s = u['src']
    loc = (f"archive: {s.get('archive')}\npages: {s.get('pages') or [s.get('page')]}\n"
           f"original_question: {s.get('original_question')}\nurl: {s.get('url')}")
    fig_lines = ''.join(f"- {u['dir']}/figures/{f}.png\n" for f in figs) or '- none\n'
    open(f"{out}/packets/{slug}.md", 'w').write(
        f"# {u['id']}\n\n## Source locator\n```\n{loc}\n```\n\n## Question figures you may open\n{fig_lines}\n"
        f"## Question (verbatim question.tex)\n```latex\n{q.strip()}\n```\n")
    rows.append(dict(id=u['id'], unit_dir=u['dir'], packet=f'packets/{slug}.md', rework=f'rework/{slug}.md'))
yaml.safe_dump(dict(course=course, population=len(units), strata=strata, seed=seed,
                    method='random.Random(seed).sample per stratum (sheet, then tripos), each ceil(20%), population sorted by id; population = ai-draft units with course or also_courses equal to the course',
                    sample=rows), open(out + '/sample.yaml', 'w'), sort_keys=False, allow_unicode=True)
print(course, 'population', len(units), strata, 'sample', len(rows))
