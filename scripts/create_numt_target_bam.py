#!/usr/bin/env python3
"""Create a compact nuclear-NUMT+chrM BAM, keeping the original intact.

Target BAM contains alignments overlapping reference NUMT loci (+padding) or
chrM. Retains read sequences and alignment flags from the original BAM for
later IGV/manual audit and validation reruns within the same regions.

NO original BAM deletion. This is not a lossless whole-genome archive. It may
omit mates/supplementary records mapped elsewhere, so full re-analysis should
use the original BAM while it is retained.
"""
import argparse
import csv
import json
from pathlib import Path
import pysam


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--loci", required=True, help="chrM_candidate_locus_validation.tsv")
    p.add_argument("--out", required=True, help="Output target BAM (different from input)")
    p.add_argument("--padding", type=int, default=1000)
    p.add_argument("--chrm", default="chrM")
    p.add_argument("--threads", type=int, default=2)
    args = p.parse_args()
    if args.padding < 100 or args.threads < 1:
        p.error("padding must be >=100 and threads >=1")
    source, dest = Path(args.bam), Path(args.out)
    if source.resolve() == dest.resolve():
        p.error("Output must not overwrite input BAM")
    if not source.is_file() or not Path(str(source)+".bai").is_file():
        p.error("Requires existing indexed BAM")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".partial")
    regions = []
    with open(args.loci, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            regions.append((r["nuclear_chrom"], int(r["nuclear_start0"]),
                            int(r["nuclear_end0"])))
    if not regions:
        p.error("No candidate NUMT loci found")
    with pysam.AlignmentFile(str(source),"rb",reference_filename=args.reference) as input_bam:
        refs = {c: input_bam.get_reference_length(c) for c in input_bam.references}
        if args.chrm not in refs:
            p.error("BAM has no chrM contig")
        segments = {args.chrm: [(0,refs[args.chrm])]}
        for chrom,start,end in regions:
            if chrom not in refs or not (0 <= start < end <= refs[chrom]):
                raise ValueError(f"Invalid locus interval: {chrom}:{start}-{end}")
            segments.setdefault(chrom,[]).append((max(0,start-args.padding),
                                                   min(refs[chrom],end+args.padding)))
        merged = []
        for chrom, rr in segments.items():
            for start,end in sorted(rr):
                if merged and merged[-1][0]==chrom and start<=merged[-1][2]:
                    merged[-1]=(chrom,merged[-1][1],max(merged[-1][2],end))
                else:
                    merged.append((chrom,start,end))
        # Unique records may occur in disjoint intervals of a chimeric read;
        # samtools view -L de-duplicates alignments across overlapping BED regions.
        # pysam loop over merged intervals does not if separate intervals touch the
        # same long alignment, so deduplicate using stable alignment coordinates.
        seen=set()
        n=0
        with pysam.AlignmentFile(str(tmp),"wb",template=input_bam) as out_bam:
            for chrom,start,end in merged:
                for read in input_bam.fetch(chrom,start,end):
                    key=(read.query_name,read.flag,read.reference_id,read.reference_start,
                         read.cigarstring,read.next_reference_id,read.next_reference_start)
                    if key in seen:
                        continue
                    seen.add(key)
                    out_bam.write(read)
                    n+=1
    # Writing reads by genomic interval means header order may be broken between
    # chromosomes. Sort explicitly before indexing.
    sorted_bam=dest.with_suffix(dest.suffix+".sorted.partial")
    try:
        pysam.sort("-@",str(args.threads),"-o",str(sorted_bam),str(tmp))
        with pysam.AlignmentFile(str(sorted_bam),"rb") as f:
            if f.check_index():
                pass
        # check_index requires index; index only after promoting the sorted BAM.
    except Exception:
        # check_index may raise because index is not yet generated; handle below.
        if not sorted_bam.exists() or sorted_bam.stat().st_size == 0:
            raise
    sorted_bam.replace(dest)
    pysam.index("-@",str(args.threads),str(dest))
    with pysam.AlignmentFile(str(dest),"rb") as bam:
        if not bam.has_index():
            raise RuntimeError("Target BAM index verification failed")
    tmp.unlink(missing_ok=True)
    metadata=dict(input=str(source),output=str(dest),padding=args.padding,
                  n_numt_loci=len(regions),n_intervals=len(merged),n_bam_records=n,
                  original_bam_preserved=True)
    dest.with_suffix(dest.suffix+".audit.json").write_text(json.dumps(metadata,indent=2)+"\n")
    print(json.dumps(metadata,indent=2))


if __name__=="__main__":
    main()
