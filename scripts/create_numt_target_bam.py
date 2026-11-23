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
    regions = []
    with open(args.loci, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            regions.append((r["nuclear_chrom"], int(r["nuclear_start0"]),
                            int(r["nuclear_end0"])))
    if not regions:
        p.error("No candidate NUMT loci found")
    with pysam.AlignmentFile(str(source),"rb",reference_filename=args.reference) as bam:
        refs = {c: bam.get_reference_length(c) for c in bam.references}
    if args.chrm not in refs:
        p.error("BAM has no chrM contig")
    segments = {args.chrm: [(0, refs[args.chrm])]}
    for chrom, start, end in regions:
        if chrom not in refs or not (0 <= start < end <= refs[chrom]):
            raise ValueError(f"Invalid locus interval: {chrom}:{start}-{end}")
        segments.setdefault(chrom, []).append(
            (max(0,start-args.padding),min(refs[chrom],end+args.padding)))
    merged=[]
    for chrom,rr in segments.items():
        for start,end in sorted(rr):
            if merged and merged[-1][0]==chrom and start<=merged[-1][2]:
                merged[-1]=(chrom,merged[-1][1],max(merged[-1][2],end))
            else:
                merged.append((chrom,start,end))
    bed=dest.with_suffix(dest.suffix+".regions.bed")
    raw=dest.with_suffix(dest.suffix+".unsorted.partial.bam")
    sorted_bam=dest.with_suffix(dest.suffix+".sorted.partial.bam")
    try:
        with bed.open("w") as handle:
            for chrom,start,end in merged:
                handle.write(f"{chrom}\\t{start}\\t{end}\\n")
        # samtools -L emits each alignment at most once even for BED overlap.
        pysam.view("-b","-L",str(bed),"-o",str(raw),str(source),catch_stdout=False)
        pysam.sort("-@",str(args.threads),"-o",str(sorted_bam),str(raw))
        with pysam.AlignmentFile(str(sorted_bam),"rb") as bam:
            n=bam.count(until_eof=True)
        sorted_bam.replace(dest)
        pysam.index("-@",str(args.threads),str(dest))
        with pysam.AlignmentFile(str(dest),"rb") as bam:
            if not bam.has_index():
                raise RuntimeError("Target BAM index verification failed")
    finally:
        bed.unlink(missing_ok=True)
        raw.unlink(missing_ok=True)
        sorted_bam.unlink(missing_ok=True)
    metadata=dict(input=str(source),output=str(dest),padding=args.padding,
                  n_numt_loci=len(regions),n_intervals=len(merged),n_bam_records=n,
                  original_bam_preserved=True)
    dest.with_suffix(dest.suffix+".audit.json").write_text(json.dumps(metadata,indent=2)+"\n")
    print(json.dumps(metadata,indent=2))


if __name__=="__main__":
    main()
