#!/usr/bin/env python3
"""Create per-run manifest; use only manually supplied local FASTQ/BAM/CRAM paths, never download ambiguous runs."""
import argparse,csv,glob
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--resolved',required=True);p.add_argument('--reference-dir',required=True);p.add_argument('--besthit-dir',required=True);p.add_argument('--local-reads',required=True,help='TSV run_accession local_path read_preset; per-run verified');p.add_argument('--out',required=True);a=p.parse_args()
reads={}
with open(a.local_reads) as f:
 for r in csv.DictReader(f,delimiter='\t'):reads[r['run_accession']]=r
out=[]
with open(a.resolved) as f:
 for r in csv.DictReader(f,delimiter='\t'):
  if r['read_type'] not in ('PACBIO','ONT') or r['run_accession'] not in reads:continue
  rp=reads[r['run_accession']];source=Path(rp['local_path']);sp=r['species'];sam=r['sample']
  refs=[p for p in glob.glob(str(Path(a.reference_dir)/f'{sp}*')) if p.endswith(('.fa','.fna','.fasta'))]
  bets=[]
  for pattern in (f'{sam}.numt_vs_chrM.besthit.tsv',f'{sam}.numt_besthit.tsv',f'{sam}.besthit.tsv'):
   bets.extend(glob.glob(str(Path(a.besthit_dir)/'**'/pattern),recursive=True))
  if not source.is_file() or not refs or not bets:continue
  preset=rp.get('read_preset','').strip() or ('map-hifi' if 'HIFI' in r['run_model'].upper() or 'REVIO' in r['run_model'].upper() or 'SEQUEL IIE' in r['run_model'].upper() else ('map-ont' if r['read_type']=='ONT' else 'map-pb'))
  out.append([sam,sp,r['run_accession'],r['read_type'],preset,str(source),sorted(refs)[0],sorted(bets)[0]])
with open(a.out,'w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['sample','species','run','platform','preset','reads','reference','besthit']);w.writerows(out)
print('ready runs:',len(out),'(missing input/reference/besthit skipped)')
