#!/usr/bin/env python3
"""Audit reference-encoded predicted NUMTs with long-read spanning support.
Coordinates: nuclear BED starts 0-based and ends exclusive; chrM BLAST 1-based inclusive.
A nuclear read spanning the reference locus is NOT by itself proof of mtDNA origin.
"""
import argparse,csv,subprocess
from pathlib import Path
import pysam
p=argparse.ArgumentParser();p.add_argument('--bam',required=True);p.add_argument('--reference',required=True);p.add_argument('--besthit',required=True);p.add_argument('--sample',required=True);p.add_argument('--run',required=True);p.add_argument('--outdir',required=True);p.add_argument('--flank',type=int,default=500);p.add_argument('--min-mapq',type=int,default=20);p.add_argument('--min-spanning',type=int,default=2);a=p.parse_args()
output=Path(a.outdir);output.mkdir(parents=True,exist_ok=True)
fa=pysam.FastaFile(a.reference);bam=pysam.AlignmentFile(a.bam,'rb',reference_filename=a.reference)
with open(a.besthit) as f: best=list(csv.DictReader(f,delimiter='\t'))
with open(output/'nuclear_support.tsv','w',newline='') as f, open(output/'reference_numt.fa','w') as fasta:
 cols=['sample','run','locus_id','nuclear_chrom','nuclear_start0','nuclear_end0','short_chrM_start','short_chrM_end','short_pident','short_strand','n_spanning','n_left_junction','n_right_junction','n_overlapping','coverage_assessable','reference_segment_status']
 w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader()
 for k,r in enumerate(best):
  try: chrom=r['nuclear_chrom'];start=int(float(r['nuclear_start']));end=int(float(r['nuclear_end']));assert chrom in fa.references and start>=0 and end>start
  except (ValueError,KeyError,AssertionError):continue
  left=set();right=set();both=set();overlap=set()
  for read in bam.fetch(chrom,max(0,start-a.flank),min(fa.get_reference_length(chrom),end+a.flank)):
   if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate or read.mapping_quality<a.min_mapq:continue
   x=read.reference_start;y=read.reference_end or x
   key=read.query_name
   if x<end and y>start:overlap.add(key)
   if x<=start-a.flank and y>=start+a.flank:left.add(key)
   if x<=end-a.flank and y>=end+a.flank:right.add(key)
   if x<=start-a.flank and y>=end+a.flank:both.add(key)
  locus=r.get('locus_id') or f'{chrom}:{start}-{end}'
  seq=fa.fetch(chrom,start,end)
  fid=f'NUMT_{k:06d}'
  fasta.write(f'>{fid} sample={a.sample} locus={locus}\n{seq}\n')
  w.writerow({'sample':a.sample,'run':a.run,'locus_id':fid,'nuclear_chrom':chrom,'nuclear_start0':start,'nuclear_end0':end,'short_chrM_start':r.get('chrM_hit_start',r.get('chrM_start','')),'short_chrM_end':r.get('chrM_hit_end',r.get('chrM_end','')),'short_pident':r.get('pident',''),'short_strand':r.get('strand',''),'n_spanning':len(both),'n_left_junction':len(left),'n_right_junction':len(right),'n_overlapping':len(overlap),'coverage_assessable':int(bool(overlap)),'reference_segment_status':'SPANNING_SUPPORTED' if len(both)>=a.min_spanning else ('PARTIAL_SUPPORT' if left and right else 'NOT_VALIDATED')})
print(output/'nuclear_support.tsv')
