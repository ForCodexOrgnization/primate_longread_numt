#!/usr/bin/env python3
"""Rerun post-alignment NUMT validators in isolated, parameter-tagged folders.

Does not touch FASTQ, BAM or original per-run validation. Run one manifest row
(or --sample/--run filter) at a time, ideally through a bounded Slurm job.
"""
import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path


def nonempty(p):
    p = Path(p)
    return p.is_file() and p.stat().st_size > 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--root", default=str(Path(__file__).resolve().parent.parent))
    ap.add_argument("--sample", required=True)
    ap.add_argument("--run", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-identity", type=float, default=.80)
    ap.add_argument("--min-locus-fraction", type=float, default=.90)
    ap.add_argument("--min-flank-fraction", type=float, default=.80)
    ap.add_argument("--max-deletion", type=int, default=20)
    ap.add_argument("--min-blast-identity", type=float, default=75)
    ap.add_argument("--min-blast-length", type=int, default=75)
    ap.add_argument("--max-reads-per-locus", type=int, default=200)
    ap.add_argument("--threads", type=int, default=2)
    args = ap.parse_args()
    if not args.tag or any(not (c.isalnum() or c in "_-") for c in args.tag):
        ap.error("tag must contain only letters/numbers/_/-")
    if args.tag.lower() in ("original", "default"):
        ap.error("Use a descriptive non-default tag to avoid confusion")
    root = Path(args.root).resolve()
    with open(args.manifest, newline="") as handle:
        rows = [r for r in csv.DictReader(handle, delimiter="\t")
                if r["sample"] == args.sample and r["run"] == args.run]
    if len(rows) != 1:
        ap.error(f"Expected exactly one manifest row, found {len(rows)}")
    r = rows[0]
    base = root / "results" / args.sample / args.run
    bam = base / "aligned.bam"
    reference = Path(r["reference"])
    besthit = Path(r["besthit"])
    for file in (bam, Path(str(bam)+".bai"), reference,
                 Path(str(reference)+".fai"), besthit):
        if not nonempty(file):
            raise SystemExit(f"Missing required input: {file}")
    if not (base/"complete.flag").is_file():
        raise SystemExit(f"Alignment run not complete: {base}")
    nuc = base / ("chrm_candidate_validation_" + args.tag)
    chrm = base / ("longread_chrM_homology_" + args.tag)
    config = dict(min_mapq=args.min_mapq, min_identity=args.min_identity,
                  min_locus_fraction=args.min_locus_fraction,
                  min_flank_fraction=args.min_flank_fraction,
                  max_deletion=args.max_deletion,
                  min_blast_identity=args.min_blast_identity,
                  min_blast_length=args.min_blast_length,
                  max_reads_per_locus=args.max_reads_per_locus,
                  reference=str(reference),besthit=str(besthit))
    provenance = base / ("revalidation_" + args.tag + ".json")
    if provenance.exists():
        old = json.loads(provenance.read_text())
        if old != config:
            raise SystemExit(f"Tag {args.tag!r} already has different parameters. Choose another tag.")
        print("Existing scenario matches:",args.tag)
        return
    if nuc.exists() or chrm.exists():
        raise SystemExit(f"Scenario dirs already exist without provenance; manually inspect before rerun: {nuc} {chrm}")
    nuc.mkdir(parents=True)
    chrm.mkdir(parents=True)
    first = [sys.executable,str(root/"scripts/validate_chrm_candidates.py"),
             "--bam",str(bam),"--reference",str(reference),"--besthit",str(besthit),
             "--outdir",str(nuc),"--min-mapq",str(args.min_mapq),
             "--min-identity",str(args.min_identity),
             "--min-locus-fraction",str(args.min_locus_fraction),
             "--min-flank-fraction",str(args.min_flank_fraction),
             "--max-deletion",str(args.max_deletion)]
    subprocess.run(first,check=True)
    second = [sys.executable,str(root/"scripts/validate_read_chrm_homology.py"),
              "--bam",str(bam),"--reference",str(reference),
              "--nuclear-validation",str(nuc),"--outdir",str(chrm),
              "--max-reads-per-locus",str(args.max_reads_per_locus),
              "--min-blast-identity",str(args.min_blast_identity),
              "--min-blast-length",str(args.min_blast_length),
              "--threads",str(args.threads)]
    subprocess.run(second,check=True)
    for fp in (nuc/"chrM_candidate_locus_validation.tsv",
               nuc/"chrM_candidate_read_evidence.tsv",
               chrm/"per_locus_chrM_homology.tsv",
               chrm/"per_read_chrM_homology.tsv"):
        if not nonempty(fp):
            raise SystemExit(f"Validation output is missing: {fp}")
    provenance.write_text(json.dumps(config,indent=2)+"\n")
    print("REVALIDATION COMPLETE:",args.sample,args.run,args.tag)


if __name__=="__main__":
    main()
