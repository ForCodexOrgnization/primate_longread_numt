#!/usr/bin/env python3
"""Prepare bounded sample-platform targeted BAM merges. Explicit --submit required."""
import argparse,csv,datetime,subprocess
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
def read(path):
    with open(path,newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest",action="append",required=True)
    p.add_argument("--out-root",default="results/sample_platform_bams")
    p.add_argument("--batch-size",type=int,default=4)
    p.add_argument("--concurrent",type=int,default=2)
    p.add_argument("--allow-partial",action="store_true")
    p.add_argument("--submit",action="store_true")
    p.add_argument("--cleanup-if-validated",action="store_true",help="After a successful full merge, delete per-run whole BAMs only if strict >90% gate passes")
    p.add_argument("--min-candidates",type=int,default=10)
    a=p.parse_args()
    if min(a.batch_size,a.concurrent,a.min_candidates)<1:p.error("counts must be >=1")
    if a.cleanup_if_validated and a.allow_partial:
        p.error("Cleanup is incompatible with --allow-partial")
    groups=defaultdict(dict)
    manifests=[str((ROOT/Path(x)).resolve()) for x in a.manifest]
    for path in manifests:
        for r in read(path):
            k=(r["sample"],r["platform"])
            if r["run"] in groups[k] and groups[k][r["run"]]!=r:
                raise ValueError(f"Conflicting manifest row: {k}/{r['run']}")
            groups[k][r["run"]]=r
    out=(ROOT/a.out_root).resolve()
    queued=[];skipped=[]
    for (sample,platform),runs in sorted(groups.items()):
        dest=out/sample/platform/"numt_target_merged.bam"
        audit=out/sample/platform/"numt_target_merged.audit.json"
        if dest.exists() or audit.exists():
            skipped.append((sample,platform,"exists"));continue
        missing=[]
        for run in runs:
            base=ROOT/"results"/sample/run
            required=[base/"complete.flag",base/"aligned.bam",base/"aligned.bam.bai",
                      base/"chrm_candidate_validation"/"chrM_candidate_locus_validation.tsv"]
            if any(not z.is_file() or (z.name!="complete.flag" and z.stat().st_size==0) for z in required):
                missing.append(run)
        if missing and not a.allow_partial:
            skipped.append((sample,platform,f"{len(missing)} runs not ready"));continue
        if len(queued)<a.batch_size:queued.append((sample,platform))
    print(f"groups={len(groups)} queued={len(queued)} skipped={len(skipped)}")
    for sample,platform in queued:print("QUEUE",sample,platform,len(groups[(sample,platform)]),"runs")
    for sample,platform,why in skipped[:20]:print("SKIP",sample,platform,why)
    if not queued:return
    if not a.submit:
        print("PLAN ONLY: add --submit after checking disk space and active Slurm jobs.")
        print("Cleanup after merge:",a.cleanup_if_validated)
        return
    folder=ROOT/"data"/"batches";folder.mkdir(parents=True,exist_ok=True)
    batch=folder/("numt_merge_"+datetime.datetime.now().strftime("%Y%m%d_%H%M%S_%f")+".tsv")
    with batch.open("w",newline="") as f:
        wr=csv.writer(f,delimiter="\t",lineterminator="\n")
        wr.writerow(["sample","platform"]);wr.writerows(queued)
    command=["sbatch",f"--array=1-{len(queued)}%{min(a.concurrent,len(queued))}",
             "--export=ALL,LR_PROJECT_ROOT="+str(ROOT)+",NUMT_MERGE_BATCH="+str(batch)+
             ",NUMT_MERGE_MANIFESTS="+",".join(manifests)+
             ",NUMT_MERGE_OUT_ROOT="+str(out)+
             ",NUMT_MERGE_ALLOW_PARTIAL="+("1" if a.allow_partial else "0")+
             ",NUMT_MERGE_CLEANUP="+("1" if a.cleanup_if_validated else "0")+
             ",NUMT_MERGE_MIN_CANDIDATES="+str(a.min_candidates),
             str(ROOT/"slurm"/"merge_sample_platform_numt.sbatch")]
    print(" ".join(command))
    subprocess.run(command,check=True,cwd=ROOT)
if __name__=="__main__":main()
