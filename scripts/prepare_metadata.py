#!/usr/bin/env python3
import argparse,csv
from pathlib import Path

def main():
 p=argparse.ArgumentParser();p.add_argument('--metadata',required=True);p.add_argument('--outdir',required=True)
 p.add_argument('--include-longread-only',action='store_true',help='Include all cohort samples with ONT or PacBio, even without Illumina; default retains existing paired-cohort filter')
 a=p.parse_args()
 out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
 with open(a.metadata,newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
 selected=[r for r in rows if ('PACBIO' in r['sequencing_platform'].upper() or 'NANOPORE' in r['sequencing_platform'].upper()) and (a.include_longread_only or 'ILLUMINA' in r['sequencing_platform'].upper())]
 cols=['input_accession','canonical_sample_id','biosample_accession','pipeline_species','sequencing_platform','sequencing_instrument_model','run_accessions','archive_mapping_n_runs']
 sample_file='all_longread_samples.tsv' if a.include_longread_only else 'matched_samples.tsv'
 run_file='all_longread_candidate_runs.tsv' if a.include_longread_only else 'candidate_runs.tsv'
 with open(out/sample_file,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows({k:r.get(k,'') for k in cols} for r in selected)
 with open(out/run_file,'w',newline='') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['sample','biosample','species','run_accession','sample_platforms','sample_instruments'])
  for r in selected:
   for accession in r['run_accessions'].split(';'):
    if accession.strip():w.writerow([r['input_accession'],r['biosample_accession'],r['pipeline_species'],accession.strip(),r['sequencing_platform'],r['sequencing_instrument_model']])
 print(f'{len(rows)} metadata rows; selected {len(selected)} samples, {len(set(x["pipeline_species"] for x in selected))} species; files: {sample_file}, {run_file}')
if __name__=='__main__':main()
