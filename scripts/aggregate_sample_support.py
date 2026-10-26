#!/usr/bin/env python3
"""Aggregate NUMT locus support by biological sample and technology.
Runs are not independent samples. Cross-run sums are NOT molecule-deduplicated.
"""
import argparse,csv,json
from collections import defaultdict
from pathlib import Path

METRICS=('n_spanning','n_left_junction','n_right_junction','n_overlapping')
NEEDED=('complete.flag','aligned.bam','aligned.bam.bai','nuclear_support.tsv','run_provenance.json')

def read_table(path):
 with open(path,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def write_table(path,rows,columns):
 path.parent.mkdir(parents=True,exist_ok=True)
 with open(path,'w',newline='') as f:
  w=csv.DictWriter(f,delimiter='\t',fieldnames=columns,extrasaction='ignore');w.writeheader();w.writerows(rows)

def main():
 p=argparse.ArgumentParser()
 p.add_argument('--manifest',default='data/remote_manifest.tsv')
 p.add_argument('--results-root',default='results')
 p.add_argument('--outdir',default='results/sample_aggregation')
 p.add_argument('--min-spanning',type=int,default=2)
 a=p.parse_args()
 runs=read_table(a.manifest);root=Path(a.results_root);out=Path(a.outdir)
 unique=set();statuses=[];details=[];expected=defaultdict(set);complete=defaultdict(set)
 for r in runs:
  sam,run=r['sample'],r['run'];key=(sam,run)
  if key in unique:raise ValueError('Duplicate sample/run: '+str(key))
  unique.add(key);expected[sam].add(run)
  d=root/sam/run
  missing=[name for name in NEEDED if not (d/name).is_file() or (name!='complete.flag' and (d/name).stat().st_size==0)]
  ok=not missing
  statuses.append(dict(sample=sam,run=run,platform=r['platform'],preset=r['preset'],status='COMPLETE' if ok else 'INCOMPLETE',missing_files=';'.join(missing)))
  if not ok:continue
  provenance=json.loads((d/'run_provenance.json').read_text())
  if provenance.get('run')!=run:raise ValueError('Provenance mismatch: '+str(key))
  complete[sam].add(run)
  for locus in read_table(d/'nuclear_support.tsv'):
   if locus['sample']!=sam or locus['run']!=run:raise ValueError('Sample/run mismatch: '+str(key))
   site=(locus['nuclear_chrom'],locus['nuclear_start0'],locus['nuclear_end0'],locus.get('short_chrM_start',''),locus.get('short_chrM_end',''))
   details.append((sam,r['platform'],run,site,locus))
 write_table(out/'run_status.tsv',statuses,['sample','run','platform','preset','status','missing_files'])
 groups=defaultdict(list)
 for sam,platform,run,site,locus in details:groups[(sam,platform,site)].append((run,locus))
 platform_rows=[]
 for (sam,platform,site),items in sorted(groups.items()):
  values={metric:[int(row.get(metric,0) or 0) for _,row in items] for metric in METRICS}
  maxes={metric:max(nums) for metric,nums in values.items()}
  status=('SPANNING_SUPPORTED' if maxes['n_spanning']>=a.min_spanning else
          'PARTIAL_SUPPORT' if maxes['n_left_junction']>0 and maxes['n_right_junction']>0 else 'NOT_VALIDATED')
  platform_rows.append(dict(sample=sam,platform=platform,nuclear_chrom=site[0],nuclear_start0=site[1],nuclear_end0=site[2],
    short_chrM_start=site[3],short_chrM_end=site[4],n_completed_runs=len(items),runs=';'.join(sorted(x[0] for x in items)),
    **{m+'_max_per_run':maxes[m] for m in METRICS},
    **{m+'_sum_runs_not_deduplicated':sum(values[m]) for m in METRICS},support_status=status))
 base=['sample','platform','nuclear_chrom','nuclear_start0','nuclear_end0','short_chrM_start','short_chrM_end','n_completed_runs','runs']
 metrics=[m+s for m in METRICS for s in ('_max_per_run','_sum_runs_not_deduplicated')]
 write_table(out/'sample_platform_locus_support.tsv',platform_rows,base+metrics+['support_status'])
 by_locus=defaultdict(list)
 for row in platform_rows:
  key=tuple(row[k] for k in ('sample','nuclear_chrom','nuclear_start0','nuclear_end0','short_chrM_start','short_chrM_end'))
  by_locus[key].append(row)
 locus_rows=[]
 for key,items in sorted(by_locus.items()):
  sam=key[0];state={r['support_status'] for r in items}
  status='SPANNING_SUPPORTED' if 'SPANNING_SUPPORTED' in state else ('PARTIAL_SUPPORT' if 'PARTIAL_SUPPORT' in state else 'NOT_VALIDATED')
  locus_rows.append(dict(zip(('sample','nuclear_chrom','nuclear_start0','nuclear_end0','short_chrM_start','short_chrM_end'),key),
   platforms=';'.join(sorted(r['platform'] for r in items)),n_completed_runs=len(complete[sam]),n_expected_runs=len(expected[sam]),
   sample_complete=int(complete[sam]==expected[sam]),max_spanning_per_run=max(r['n_spanning_max_per_run'] for r in items),
   sum_spanning_across_runs_not_deduplicated=sum(r['n_spanning_sum_runs_not_deduplicated'] for r in items),support_status=status))
 locuscols=['sample','nuclear_chrom','nuclear_start0','nuclear_end0','short_chrM_start','short_chrM_end','platforms','n_completed_runs','n_expected_runs','sample_complete','max_spanning_per_run','sum_spanning_across_runs_not_deduplicated','support_status']
 write_table(out/'sample_locus_support.tsv',locus_rows,locuscols)
 summaries=[]
 for sam in sorted(expected):
  loci=[r for r in locus_rows if r['sample']==sam]
  summaries.append(dict(sample=sam,n_expected_runs=len(expected[sam]),n_completed_runs=len(complete[sam]),
   n_pending_runs=len(expected[sam]-complete[sam]),sample_complete=int(expected[sam]==complete[sam]),
   n_loci_with_results=len(loci),n_spanning_supported=sum(r['support_status']=='SPANNING_SUPPORTED' for r in loci),
   n_partial_support=sum(r['support_status']=='PARTIAL_SUPPORT' for r in loci),
   n_not_validated=sum(r['support_status']=='NOT_VALIDATED' for r in loci)))
 write_table(out/'sample_summary.tsv',summaries,list(summaries[0]) if summaries else ['sample'])
 print(json.dumps(dict(samples=len(expected),completed_runs=sum(len(v) for v in complete.values()),sample_loci=len(locus_rows)),indent=2))
if __name__=='__main__':main()
