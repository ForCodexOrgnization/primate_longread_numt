#!/usr/bin/env python3
"""Create run-level tasks from ENA-verified metadata; fail closed on ambiguity."""
import argparse,csv,glob,re
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--resolved',required=True);p.add_argument('--reference-dir',required=True)
p.add_argument('--besthit-dir',required=True);p.add_argument('--out',default='data/remote_manifest.tsv')
p.add_argument('--presets');a=p.parse_args()
presets={}
if a.presets:
 with open(a.presets) as f: presets={r['run_accession']:r['read_preset'] for r in csv.DictReader(f,delimiter='\t')}
ready=[];skipped=[]
with open(a.resolved) as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r.get('read_type') not in ('PACBIO','ONT'): continue
  sample,run,species=r['sample'],r['run_accession'],r['species']
  refs=sorted(x for x in glob.glob(str(Path(a.reference_dir)/(species+'*'))) if x.endswith(('.fa','.fna','.fasta')) and Path(x).is_file())
  hits=set()
  for suffix in ('.numt_vs_chrM.besthit.tsv','.numt_besthit.tsv','.besthit.tsv'):
   hits.update(glob.glob(str(Path(a.besthit_dir)/'**'/(sample+suffix)),recursive=True))
  hits=sorted(x for x in hits if Path(x).is_file())
  urls=[x.strip() for x in r.get('fastq_ftp','').split(';') if x.strip()]
  valid=bool(urls) and all(re.fullmatch(r'(?:ftp\.sra\.ebi\.ac\.uk/|https://ftp\.sra\.ebi\.ac\.uk/)[A-Za-z0-9_./+%-]+',u) and u.lower().endswith(('.fastq.gz','.fq.gz','.fastq','.fq')) for u in urls)
  preset=presets.get(run,'map-ont' if r['read_type']=='ONT' else '')
  reasons=[]
  if len(refs)!=1: reasons.append('reference_missing_or_ambiguous')
  if len(hits)!=1: reasons.append('besthit_missing_or_ambiguous')
  if not valid: reasons.append('fastq_ftp_missing_or_invalid')
  if preset not in ('map-pb','map-hifi','map-ont','lr:hq'): reasons.append('preset_requires_review')
  if reasons: skipped.append([run,sample,';'.join(reasons)]);continue
  ready.append([sample,species,run,r['read_type'],preset,';'.join(urls),refs[0],hits[0]])
out=Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
with out.open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['sample','species','run','platform','preset','fastq_ftp','reference','besthit']);w.writerows(ready)
with out.with_suffix('.skipped.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['run','sample','reason']);w.writerows(skipped)
print('Ready runs:',len(ready),'Skipped:',len(skipped))
