#!/usr/bin/env python3
"""Restartable, bounded, Slurm-aware long-read cohort scheduler.

One controller process owns the state lock; call --once to preview or with
--submit to dispatch. --watch periodically polls/dispatches from a compute
node. Only aria2 checkpoint exit 75, TIMEOUT, PREEMPTED, NODE_FAIL,
BOOT_FAIL trigger automatic resubmission, at most --max-attempts attempts.
All other errors block for inspection. Missing sacct accounting BLOCKS (never
guesses job failure). Checks active Slurm queues before every submission.
"""
import argparse
import csv
import fcntl
import json
import os
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
FIELDS=["sample","species","run","platform","preset","fastq_ftp","reference","besthit"]
GOOD_STATES={"COMPLETED"}
ACTIVE={"PENDING","RUNNING","CONFIGURING","COMPLETING","SUSPENDED","RESIZING","REQUEUED","REQUEUE_FED","STAGE_OUT"}
RETRY_STATES={"TIMEOUT","NODE_FAIL","PREEMPTED","BOOT_FAIL"}
MIN_FREE_MULTIPLIER=1.0


def stamp():
    return datetime.now(timezone.utc).isoformat()


def read_tsv(path):
    with open(path,newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))


def isfile(path):
    p=Path(path)
    return p.is_file() and p.stat().st_size>0


def output_complete(run,stage,results):
    base=results/run["sample"]/run["run"]
    if stage=="align":
        return ((base/"complete.flag").is_file() and all(isfile(base/x) for x in
                ("aligned.bam","aligned.bam.bai","nuclear_support.tsv","reference_numt.fa","run_provenance.json")))
    if stage=="validate":
        return all(isfile(base/x) for x in
                   ("chrm_candidate_validation/chrM_candidate_locus_validation.tsv",
                    "chrm_candidate_validation/chrM_candidate_read_evidence.tsv",
                    "longread_chrM_homology/per_locus_chrM_homology.tsv",
                    "longread_chrM_homology/per_read_chrM_homology.tsv"))
    raise ValueError(stage)


def cmd(args):
    p=subprocess.run([str(x) for x in args],capture_output=True,text=True)
    if p.returncode:
        raise RuntimeError(f"Command failed: {args}: {p.stderr.strip() or p.stdout.strip()}")
    return p.stdout.strip()


def query_slurm(jobid):
    """Return running / terminal / unknown, with state and exit code."""
    queue=cmd(["squeue","-h","-j",str(jobid),"-o","%T"])
    if queue.strip():
        return ("active",queue.splitlines()[0].strip().upper(),None)
    # Allow for squeue-to-sacct propagation. Missing data is NOT failure.
    accounting=cmd(["sacct","-n","-P","-j",str(jobid),
                    "-o","JobIDRaw,State,ExitCode"])
    for line in accounting.splitlines():
        fields=line.split("|")
        if len(fields)<3 or fields[0]!=str(jobid):
            continue
        state=fields[1].split(" ",1)[0].split("+",1)[0].upper()
        exitcode=fields[2].split(":",1)[0]
        try: exitcode=int(exitcode)
        except ValueError: exitcode=None
        return ("terminal",state,exitcode)
    return ("unknown","NO_ACCOUNTING",None)


def group_key(row):
    return (row["sample"],row["platform"])


def load_manifest(paths):
    records={}
    for path in paths:
        for row in read_tsv(path):
            if not all(row.get(f) for f in FIELDS):
                raise ValueError(f"Incomplete manifest entry in {path}: {row}")
            if row["platform"] not in ("ONT","PACBIO"):
                raise ValueError(f"Bad platform: {row['platform']}")
            if row["preset"] not in ("map-ont","map-pb","map-hifi","lr:hq"):
                raise ValueError(f"Unreviewed preset for {row['run']}")
            if not re.fullmatch("[A-Za-z0-9._-]+",row["run"]):
                raise ValueError("Invalid run accession")
            key=row["run"]
            if key in records and records[key]!=row:
                raise ValueError(f"Conflicting run row {key}")
            records[key]=row
    if not records:
        raise SystemExit("Empty cohort manifest")
    return records


def atomic_write(path, value):
    temp=path.with_suffix(".tmp")
    with temp.open("w") as f:
        json.dump(value,f,indent=2,sort_keys=True)
        f.write("\n")
        f.flush()
        os.fsync(f.fileno())
    os.replace(temp,path)


def status_key(stage,run):
    return stage+":"+run


def submit(stage,row,root,results,state_dir):
    """One-run manifest per job, stable job name for human inspection."""
    runid=row["run"]
    manifest=state_dir/"jobs"/(stage+"_"+runid+".tsv")
    manifest.parent.mkdir(parents=True,exist_ok=True)
    with manifest.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,delimiter="\t",lineterminator="\n")
        w.writeheader()
        w.writerow({x:row[x] for x in FIELDS})
    jobname=("ln_al_" if stage=="align" else "ln_va_")+runid
    exports=["ALL","LR_PROJECT_ROOT="+str(root)]
    if stage=="align":
        exports+=["BATCH_MANIFEST="+str(manifest),"LR_OUT_ROOT="+str(results)]
        script=root/"slurm"/"download_align_validate.sbatch"
    else:
        exports+=["VALIDATION_MANIFEST="+str(manifest)]
        script=root/"slurm"/"validate_pilot_array.sbatch"
    output=cmd(["sbatch","--parsable","--job-name="+jobname,"--array=1-1",
                "--export="+",".join(exports),str(script)])
    jobid=output.split(";",1)[0]
    if not jobid.isdigit():
        raise RuntimeError("Unexpected sbatch job id: "+output)
    return jobid+"_1"


def check_queue_for_run(run,stage):
    """Prevent duplicate tasks if job was submitted outside controller."""
    jobname=("ln_al_" if stage=="align" else "ln_va_")+run
    data=cmd(["squeue","-h","-u",os.environ["USER"],"-o","%j|%i"])
    for line in data.splitlines():
        name,_,jobid=line.partition("|")
        if name==jobname:
            return jobid
    return None


def gib_free(path):
    s=os.statvfs(path)
    return s.f_bavail*s.f_frsize/(1024**3)


def run_once(a,do_submit,root,manifest,state_dir):
    results=(root/"results").resolve()
    staging=(root/"data"/"longread_download_stage").resolve()
    state_dir.mkdir(parents=True,exist_ok=True)
    (state_dir/"jobs").mkdir(exist_ok=True)
    state_path=state_dir/"state.json"
    with (state_dir/"controller.lock").open("a") as lf:
        fcntl.flock(lf,fcntl.LOCK_EX|fcntl.LOCK_NB)
        state=json.loads(state_path.read_text()) if state_path.exists() else {"jobs":{},"created_utc":stamp()}
        jobs=state.setdefault("jobs",{})
        state["last_tick_utc"]=stamp()
        state["manifest_runs"]=len(manifest)
        # After sample-platform BAM cleanup the run-level full BAM is removed,
        # but validated per-read TSVs still certify completed processing.
        groups=defaultdict(list)
        for row in manifest.values():
            groups[group_key(row)].append(row)
        order=sorted(groups,key=lambda k:(0 if k[1]=="ONT" else 1,len(groups[k]),k))
        tasks=[]
        for g in order:
            tasks+=sorted(groups[g],key=lambda r:r["run"])
        states={}
        for row in tasks:
            accession=row["run"]
            for stage in ("align","validate"):
                key=status_key(stage,accession)
                info=jobs.setdefault(key,{"attempts":0,"status":"PENDING"})
                if "status" not in info: info["status"]="PENDING"
                if info["status"]=="SUBMITTING" and not info.get("jobid"):
                    found=check_queue_for_run(accession,stage)
                    if found:
                        info.update(jobid=found,status="ACTIVE",attempts=info["attempts"]+1)
                    else:
                        info["status"]="WAIT_SUBMISSION_RECONCILIATION"
                # A completed validation is final even when whole BAM
                # has already been archived and cleaned.
                done=output_complete(row,stage,results)
                if stage=="align" and not done and output_complete(row,"validate",results):
                    done=True
                if done:
                    info["status"]="DONE"
                    states[key]="DONE"
                    continue
                if info.get("status")=="DONE":
                    # Previously completed outputs disappeared unexpectedly:
                    # fail closed rather than silently remapping.
                    info["status"]="BLOCKED_OUTPUT_MISSING"
                    states[key]=info["status"]
                    continue
                if info.get("jobid"):
                    kind,slurm,exitcode=query_slurm(info["jobid"])
                    info.update(last_slurm=slurm,last_exit_code=exitcode,
                                last_observed_utc=stamp())
                    if kind=="active":
                        info["status"]="ACTIVE"
                    elif kind=="unknown":
                        info["status"]="WAIT_ACCOUNTING"
                    else:
                        info.pop("jobid")
                        if slurm=="COMPLETED":
                            info["status"]="BLOCKED_OUTPUT_MISSING"
                        elif slurm in RETRY_STATES or (stage=="align" and exitcode==75):
                            info["status"]="RETRYABLE"
                            info["last_error"]=slurm+":"+str(exitcode)
                            info["next_retry_time"]=time.time()+a.retry_delay_minutes*60
                        else:
                            info["status"]="BLOCKED_FAILURE"
                            info["last_error"]=slurm+":"+str(exitcode)
                states[key]=info.get("status","PENDING")
        # Guard against all active pipeline work, not only our own jobs.
        counts={stage:sum(states.get(status_key(stage,r["run"])) in
                          ("ACTIVE","WAIT_ACCOUNTING") for r in tasks)
                for stage in ("align","validate")}
        selected=[]
        free=min(gib_free(results),gib_free(staging)) if staging.exists() else gib_free(results)
        disk_ok=free>=a.min_free_gib
        for row in tasks:
            for stage in ("validate","align"):
                key=status_key(stage,row["run"])
                info=jobs[key]
                if info["status"] not in ("PENDING","RETRYABLE"):
                    continue
                if stage=="validate" and not output_complete(row,"align",results):
                    continue
                if stage=="align" and (not disk_ok or counts[stage]>=a.max_align):
                    continue
                if stage=="validate" and counts[stage]>=a.max_validate:
                    continue
                if info["attempts"]>=a.max_attempts:
                    info["status"]="MAX_ATTEMPTS"
                    continue
                if info.get("next_retry_time",0)>time.time():
                    continue
                external=check_queue_for_run(row["run"],stage)
                if external:
                    info.update(jobid=external,
                                status="ACTIVE",external_discovered=True)
                    counts[stage]+=1
                    continue
                selected.append((stage,row))
                counts[stage]+=1
                info["status"]="PLANNED" if not do_submit else "SUBMITTING"
                if do_submit:
                    # Persist the intent before external submission. If a
                    # controller crashes during sbatch, the next cycle checks
                    # the stable Slurm job name before any resubmission.
                    atomic_write(state_path,state)
                    jobid=submit(stage,row,root,results,state_dir)
                    info.update(jobid=jobid,attempts=info["attempts"]+1,
                                submitted_utc=stamp(),status="ACTIVE")
                    info.pop("next_retry_time",None)
                    atomic_write(state_path,state)
            # Prioritize validation of already-aligned runs in group order.
        for info in jobs.values():
            if info.get("status")=="PLANNED":
                info["status"]="PENDING"
            if info.get("status")=="SUBMITTING":
                info["status"]="WAIT_SUBMISSION_RECONCILIATION"
        atomic_write(state_path,state)
        totals=defaultdict(int)
        for row in tasks:
            for stage in ("align","validate"):
                totals[jobs[status_key(stage,row["run"])]["status"]]+=1
        print(f"COHORT {len(tasks)} runs; disk free >= {free:.1f} GiB; min {a.min_free_gib:.1f} GiB",flush=True)
        print("JOB STATES",dict(sorted(totals.items())),flush=True)
        for stage,row in selected:
            print(("SUBMITTED" if do_submit else "WOULD_SUBMIT"),stage,row["sample"],row["run"],flush=True)
        all_done=all(jobs[status_key("validate",r["run"])]["status"]=="DONE" for r in tasks)
        blocked=sum(info["status"].startswith("BLOCKED") or info["status"]=="MAX_ATTEMPTS" for info in jobs.values())
        if not disk_ok:print("PAUSED: disk guard; no new alignments",flush=True)
        if blocked:print("ATTENTION: blocked jobs:",blocked,flush=True)
        return all_done


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest",required=True)
    ap.add_argument("--state-dir",default="data/cohort_manager")
    ap.add_argument("--min-free-gib",type=float,default=500)
    ap.add_argument("--max-align",type=int,default=1)
    ap.add_argument("--max-validate",type=int,default=1)
    ap.add_argument("--max-attempts",type=int,default=5)
    ap.add_argument("--retry-delay-minutes",type=int,default=15)
    ap.add_argument("--submit",action="store_true",help="Allow Slurm submissions; default preview")
    ap.add_argument("--watch",action="store_true",help="Poll repeatedly; run via Slurm controller job")
    ap.add_argument("--poll-seconds",type=int,default=180)
    ap.add_argument("--watch-hours",type=float,default=18)
    a=ap.parse_args()
    if not all(x>0 for x in (a.min_free_gib,a.max_align,a.max_validate,
                              a.max_attempts,a.retry_delay_minutes,
                              a.poll_seconds,a.watch_hours)):
        ap.error("All resource and timing parameters must be >0")
    if a.watch and not a.submit:
        ap.error("--watch requires --submit")
    root=ROOT
    manifest=load_manifest([(root/a.manifest).resolve()])
    state_dir=(root/a.state_dir).resolve()
    if not (root/"config/bouchet.env").is_file():
        ap.error("Run from full Bouchet checkout containing config/bouchet.env")
    stop=state_dir/"STOP"
    if stop.exists():
        print(f"STOP requested: {stop}; no new submissions")
        return
    start=time.monotonic()
    while True:
        try:
            finished=run_once(a,a.submit,root,manifest,state_dir)
        except BlockingIOError:
            print("A cohort controller already holds the lock. Exiting without submissions.")
            return
        if not a.watch or finished:
            break
        if stop.exists():
            print("STOP detected; exiting controller")
            return
        if time.monotonic()-start+a.poll_seconds>a.watch_hours*3600:
            print("WATCH_WINDOW_ENDED; restart controller or schedule next window")
            break
        time.sleep(a.poll_seconds)


if __name__=="__main__":
    main()
