#!/usr/bin/env python3
"""Select only explicitly approved pilot runs, ordered by priority.

No substitution: missing reference, besthit, FASTQ or preset remains skipped.
"""
import argparse,csv
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--resolved',default='data/resolved_runs.tsv')
p.add_argument('--selection',default='data/pilot_run_selection.tsv')
p.add_argument('--out',default='data/pilot_resolved_runs.tsv')
a=p.parse_args()
with open(a.resolved,newline='') as f: src=list(csv.DictReader(f,delimiter='\t'))
byrun={r['run_accession']:r for r in src}
with open(a.selection,newline='') as f: chosen=list(csv.DictReader(f,delimiter='\t'))
if len({r['run_accession'] for r in chosen})!=len(chosen):raise SystemExit('Duplicate run in selection')
out=[];problems=[]
for row in sorted(chosen,key=lambda x:int(x['priority'])):
 run=row['run_accession'];r=byrun.get(run)
 if r is None:problems.append((run,'not present in resolved runs'));continue
 if r['sample']!=row['sample'] or r['species']!=row['species']:
  problems.append((run,'sample/species mismatch'));continue
 if r['read_type'] not in ('ONT','PACBIO') or not r.get('fastq_ftp'):
  problems.append((run,'not verified long read or missing FASTQ'));continue
 out.append(r)
if problems:
 for x in problems:print('ERROR:',*x)
 raise SystemExit('Selection validation failed')
path=Path(a.out);path.parent.mkdir(parents=True,exist_ok=True)
with path.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(src[0]),delimiter='\t');w.writeheader();w.writerows(out)
print('Selected',len(out),'runs;',len({r['sample'] for r in out}),'samples;',len({r['species'] for r in out}),'species')
