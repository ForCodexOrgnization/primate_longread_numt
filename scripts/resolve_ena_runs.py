#!/usr/bin/env python3
"""Verify per-run sequencing platform with ENA Portal API; fail closed for unknown platforms."""
import argparse,csv,json,time,urllib.parse,urllib.request
from collections import defaultdict
from pathlib import Path

def get_run(run,retries=3):
 url='https://www.ebi.ac.uk/ena/portal/api/filereport?'+urllib.parse.urlencode({'accession':run,'result':'read_run','fields':'run_accession,sample_accession,instrument_platform,instrument_model,fastq_ftp,fastq_bytes,submitted_ftp','format':'json'})
 for i in range(retries):
  try:
   req=urllib.request.Request(url,headers={'User-Agent':'primate-longread-numt/0.1'})
   with urllib.request.urlopen(req,timeout=45) as response: data=json.load(response)
   return data[0] if data else None
  except Exception:
   if i==retries-1:raise
   time.sleep(2**i)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--candidate-runs',required=True);ap.add_argument('--out',required=True);ap.add_argument('--cache',default='');a=ap.parse_args()
 with open(a.candidate_runs) as f: rows=list(csv.DictReader(f,delimiter='\t'))
 cachepath=Path(a.cache or (a.out+'.cache.json'));cache=json.loads(cachepath.read_text()) if cachepath.exists() else {}
 output=[]
 for j,row in enumerate(rows,1):
  run=row['run_accession']
  if run not in cache:
   try:cache[run]=get_run(run)
   except Exception as e:cache[run]={'_error':str(e)}
   if j%10==0:cachepath.write_text(json.dumps(cache,indent=2))
  info=cache[run] or {}
  platform=info.get('instrument_platform','').upper()
  kind='PACBIO' if 'PACBIO' in platform else ('ONT' if 'OXFORD_NANOPORE' in platform else 'OTHER')
  output.append({**{k:row.get(k,'') for k in ['sample','biosample','species','run_accession']},'run_platform':platform,'run_model':info.get('instrument_model',''),'read_type':kind,'fastq_ftp':info.get('fastq_ftp',''),'submitted_ftp':info.get('submitted_ftp',''),'ena_sample_accession':info.get('sample_accession',''),'lookup_error':info.get('_error','')})
 cachepath.write_text(json.dumps(cache,indent=2))
 with open(a.out,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(output[0]) if output else ['sample'],delimiter='\t');w.writeheader();w.writerows(output)
 print('verified long-read runs:',sum(x['read_type'] in ('PACBIO','ONT') for x in output),'of',len(output))
if __name__=='__main__':main()
