#!/usr/bin/env python3
"""Generate a conservative, editable run-to-minimap2-preset template.
Instrument model alone cannot distinguish PacBio CLR from HiFi.
Only ONT is assigned automatically; PacBio requires verified per-run read type.
"""
import argparse,csv
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--resolved',default='data/resolved_runs.tsv')
p.add_argument('--out',default='data/presets.tsv')
a=p.parse_args()
rows=[]
with open(a.resolved,newline='') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r.get('read_type') not in ('PACBIO','ONT'):continue
  rows.append({'run_accession':r['run_accession'],'read_preset':'map-ont' if r['read_type']=='ONT' else 'REVIEW','read_type':r['read_type'],'run_model':r.get('run_model',''),'notes':'Inspect whether CCS/HiFi (map-hifi) or CLR/subreads (map-pb)' if r['read_type']=='PACBIO' else 'Default ONT preset'})
out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['run_accession','read_preset','read_type','run_model','notes'],delimiter='\t');w.writeheader();w.writerows(rows)
print('Wrote:',out,'total:',len(rows),'requires review:',sum(r['read_preset']=='REVIEW' for r in rows))
