#!/usr/bin/env python3
"""Fail-closed, optional post-merge full-BAM cleanup for one sample/platform.

Policy: strict BOTH-supported BLAST_QC_PASS loci / BLAST_QC_PASS loci > 0.90.
Computes results freshly from per-run read evidence, never trusts stale
summary tables. Requires every manifest run validated, identical candidate
locus IDs/coordinates, and a verified complete merged targeted BAM with
run-group provenance. Default is DRY RUN. --execute permanently deletes only
the approved aligned.bam and aligned.bam.bai inputs after saving an audit.
"""
import argparse
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


def read(path):
    with open(path,newline="") as handle:
        return list(csv.DictReader(handle,delimiter="\t"))


def nonempty(path):
    p=Path(path)
    return p.is_file() and p.stat().st_size>0


def run(cmd):
    return subprocess.run([str(x) for x in cmd],check=True,
                          capture_output=True,text=True).stdout


def now():
    return datetime.now(timezone.utc).isoformat()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest",action="append",required=True)
    ap.add_argument("--sample",required=True)
    ap.add_argument("--platform",choices=("ONT","PACBIO"),required=True)
    ap.add_argument("--results-root",default="results")
    ap.add_argument("--out-root",default="results/sample_platform_bams")
    ap.add_argument("--min-rate",type=float,default=0.90)
    ap.add_argument("--min-candidates",type=int,default=10,
                    help="Fail closed on very small denominators; default 10")
    ap.add_argument("--execute",action="store_true",help="Permanently delete approved whole-genome BAM/BAI")
    args=ap.parse_args()
    if not (0<args.min_rate<1) or args.min_candidates<1:
        ap.error("Invalid rate or denominator requirement")
    if not re.fullmatch(r"[A-Za-z0-9._-]+",args.sample):
        ap.error("Invalid sample")
    root=Path(args.results_root).resolve()
    directory=Path(args.out_root).resolve()/args.sample/args.platform
    bam=directory/"numt_target_merged.bam"
    bai=directory/"numt_target_merged.bam.bai"
    audit=directory/"numt_target_merged.audit.json"
    if not all(nonempty(p) for p in (bam,bai,audit)):
        raise SystemExit("SKIP: missing merged BAM, index, or merge audit")
    merged=json.loads(audit.read_text())
    if (merged.get("sample")!=args.sample or merged.get("platform")!=args.platform
        or merged.get("partial") is not False
        or merged.get("n_manifest_runs")!=merged.get("n_included_runs")):
        raise SystemExit("SKIP: merged audit does not prove a complete group")

    runs={}
    for path in args.manifest:
        for row in read(path):
            if row["sample"]!=args.sample or row["platform"]!=args.platform:
                continue
            accession=row["run"]
            if accession in runs and runs[accession]!=row:
                raise SystemExit(f"Conflicting manifest entry for {accession}")
            runs[accession]=row
    if not runs or set(runs)!=set(merged.get("included_runs",[])):
        raise SystemExit("SKIP: merged BAM source runs differ from full manifest")
    if len({r["reference"] for r in runs.values()})!=1 or len({r["besthit"] for r in runs.values()})!=1:
        raise SystemExit("SKIP: reference or NUMT candidate source mismatch")
    if (merged.get("reference")!=next(iter(runs.values()))["reference"]
        or merged.get("besthit")!=next(iter(runs.values()))["besthit"]):
        raise SystemExit("SKIP: merge audit reference/besthit mismatch")

    run_inputs={}
    candidates=None
    signatures=None
    for accession in sorted(runs):
        base=root/args.sample/accession
        nbase=base/"chrm_candidate_validation"
        cbase=base/"longread_chrM_homology"
        files=[base/"aligned.bam",base/"aligned.bam.bai",
               nbase/"chrM_candidate_locus_validation.tsv",
               nbase/"chrM_candidate_read_evidence.tsv",
               cbase/"per_locus_chrM_homology.tsv",
               cbase/"per_read_chrM_homology.tsv"]
        if not (base/"complete.flag").is_file() or not all(nonempty(f) for f in files):
            raise SystemExit(f"SKIP: incomplete alignment/validation for {accession}")
        loci=read(files[2])
        if not loci:
            raise SystemExit(f"SKIP: no candidate loci for {accession}")
        current={r["locus_id"]:(r["nuclear_chrom"],r["nuclear_start0"],
                 r["nuclear_end0"],r["chrM_hit"],r["chrM_start"],
                 r["chrM_end"],r["group"]) for r in loci}
        if len(current)!=len(loci):
            raise SystemExit(f"SKIP: duplicate locus IDs for {accession}")
        if signatures is None:
            signatures=current
        elif current!=signatures:
            raise SystemExit(f"SKIP: non-identical NUMT candidate lists across runs ({accession})")
        run_inputs[accession]=files

    # Merged read group IDs must correspond to EVERY source run.
    run(["samtools","quickcheck",bam])
    run(["samtools","idxstats",bam])
    header=run(["samtools","view","-H",bam])
    rg=set()
    for line in header.splitlines():
        if line.startswith("@RG\t"):
            fields=dict(x.split(":",1) for x in line.split("\t")[1:])
            if fields.get("SM")!=args.sample or fields.get("PL")!=args.platform:
                raise SystemExit("SKIP: merged BAM read-group sample/platform mismatch")
            rg.add(fields.get("ID"))
    if rg!=set(runs):
        raise SystemExit("SKIP: merged read groups do not match every run")

    # Recompute strict sample-platform counts from the current inputs, in a
    # private directory, so no stale TSV or older parameter scenario can pass.
    script=Path(__file__).resolve().parent/"aggregate_sample_level_validation.py"
    if not script.is_file():
        raise SystemExit("Missing sample-level aggregator script")
    with tempfile.TemporaryDirectory(prefix=".cleanup_review_",dir=directory) as td:
        work=Path(td)
        manifest=work/"group.tsv"
        cols=list(next(iter(runs.values())).keys())
        with manifest.open("w",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=cols,delimiter="\t",lineterminator="\n")
            writer.writeheader()
            writer.writerows(runs[k] for k in sorted(runs))
        command=[sys.executable,str(script),"--manifest",str(manifest),
                 "--results-root",str(root),"--outdir",str(work/"validation"),
                 "--require-complete"]
        run(command)
        values=read(work/"validation"/"sample_platform_locus_validation.tsv")
        wanted=[r for r in values if r["sample"]==args.sample and r["platform"]==args.platform]
        if len(wanted)!=len(signatures) or {r["locus_id"] for r in wanted}!=set(signatures):
            raise SystemExit("SKIP: re-aggregated loci do not match input candidates")
        passing=[r for r in wanted if r["group"]=="BLAST_QC_PASS"]
        numerator=sum(r["final_status"]=="VALIDATED" for r in passing)
        denominator=len(passing)
        rate=numerator/denominator if denominator else 0
        decision=denominator>=args.min_candidates and rate>args.min_rate
        assessment=dict(sample=args.sample,platform=args.platform,
                        assessed_utc=now(),metric="strict_both_supported_over_BLAST_QC_PASS",
                        numerator=numerator,denominator=denominator,rate=rate,
                        min_rate_exclusive=args.min_rate,min_candidates=args.min_candidates,
                        approved=decision,manifest_runs=sorted(runs),
                        merged_bam=str(bam),merged_bam_bytes=bam.stat().st_size,
                        source_bams=[str(run_inputs[k][0]) for k in sorted(runs)],
                        policy="No deletion if rate<=threshold, denominator insufficient, incomplete evidence or merged BAM not validated")
        print(json.dumps(assessment,indent=2))
        if not args.execute:
            print("DRY RUN: no files deleted; add --execute only after review")
            return
        if not decision:
            print("SKIP: validation threshold not met; no files deleted")
            return
        # Write durable decision before deletion; use exclusive create to
        # block repeated or concurrently overlapping cleanup attempts.
        record=directory/"whole_bam_cleanup_audit.json"
        fd=os.open(record,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
        with os.fdopen(fd,"w") as f:
            assessment["deletion_started_utc"]=now()
            assessment["deleted"]=[]
            json.dump(assessment,f,indent=2)
            f.write("\n")
            f.flush()
            os.fsync(f.fileno())
        deleted=[]
        for accession in sorted(runs):
            for file in run_inputs[accession][:2]:
                if not file.is_file() or file.is_symlink():
                    raise RuntimeError("Refusing missing/symlink BAM path during deletion: "+str(file))
                file.unlink()
                deleted.append(str(file))
        assessment["deleted"]=deleted
        assessment["deletion_finished_utc"]=now()
        record.write_text(json.dumps(assessment,indent=2)+"\n")
        print("DONE: removed",len(deleted),"run-level BAM/BAI files; all evidence and targeted BAM kept")


if __name__=="__main__":
    main()
