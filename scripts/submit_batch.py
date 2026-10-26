#!/usr/bin/env python3
"""Submit a bounded Slurm batch; do not overlap batches without run reservations."""
import argparse,csv,datetime,subprocess
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--manifest',default='data/remote_manifest.tsv');p.add_argument('--unit',choices=['run','sample'],default='run')
p.add_argument('--batch-size',type=int,default=5);p.add_argument('--concurrent',type=int,default=2)
p.add_argument('--max-runs',type=int,default=30);p.add_argument('--out-root',default='results');p.add_argument('--dry-run',action='store_true');a=p.parse_args()
if min(a.batch_size,a.concurrent,a.max_runs)<1:p.error('counts must be >=1')
with open(a.manifest) as f:rows=list(csv.DictReader(f,delimiter='\t'))
chosen=[];seen=set()
for r in rows:
 out=Path(a.out_root)/r['sample']/r['run']
 if (out/'complete.flag').is_file() and all((out/name).is_file() and (out/name).stat().st_size>0 for name in ('aligned.bam','aligned.bam.bai','nuclear_support.tsv','reference_numt.fa','run_provenance.json')):continue
 if a.unit=='sample' and r['sample'] not in seen and len(seen)>=a.batch_size:continue
 if a.unit=='run' and len(chosen)>=a.batch_size:break
 if len(chosen)>=a.max_runs:break
 chosen.append(r);seen.add(r['sample'])
if not chosen:print('No unfinished runs selected');raise SystemExit(0)
folder=Path('data/batches');folder.mkdir(parents=True,exist_ok=True)
manifest=folder/('batch_'+datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f')+'.tsv')
with manifest.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(chosen[0]),delimiter='\t');w.writeheader();w.writerows(chosen)
cmd=['sbatch',f'--array=1-{len(chosen)}%{min(a.concurrent,len(chosen))}','--export=ALL,BATCH_MANIFEST='+str(manifest.resolve())+',LR_OUT_ROOT='+str(Path(a.out_root).resolve()),'slurm/download_align_validate.sbatch']
print('Selected',len(chosen),'runs from',len(seen),'samples; manifest:',manifest)
print(' '.join(cmd))
if not a.dry_run:subprocess.run(cmd,check=True)
