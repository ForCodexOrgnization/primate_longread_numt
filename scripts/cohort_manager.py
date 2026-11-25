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
import hashlib
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
    # A just-finished array task may no longer exist in squeue. On
    # Bouchet, squeue -j ID_1 can exit nonzero with "Invalid job id
    # specified" instead of returning an empty listing. This is expected
    # during the queue -> accounting transition, not a controller failure.
    queue_call=subprocess.run(["squeue","-h","-j",str(jobid),"-o","%T"],
                              capture_output=True,text=True)
    if queue_call.returncode:
        message=(queue_call.stderr or "")+" "+(queue_call.stdout or "")
        if not re.search(r"Invalid job id specified",message,re.IGNORECASE):
            raise RuntimeError(f"squeue failed for {jobid}: {message.strip()}")
        queue=""
    else:
        queue=queue_call.stdout.strip()
    if queue:
        return ("active",queue.splitlines()[0].strip().upper(),None)
    # Allow for squeue-to-sacct propagation. Missing data is NOT failure.
    # Query the parent job; its accounting includes the task ID.
    parent=str(jobid).split("_",1)[0]
    accounting=cmd(["sacct","-n","-P","-j",parent,
                    "-o","JobIDRaw,State,ExitCode"])
    entries={}
    for line in accounting.splitlines():
        fields=line.split("|")
        if len(fields)<3:
            continue
        raw_id=fields[0].strip()
        state=fields[1].split(" ",1)[0].split("+",1)[0].upper()
        exitcode=fields[2].split(":",1)[0]
        try: exitcode=int(exitcode)
        except ValueError: exitcode=None
        entries[raw_id]=(state,exitcode)
    if str(jobid) in entries:
        return ("terminal",*entries[str(jobid)])
    # This controller submits exactly ONE task per array (1-1).
    # Some Bouchet sacct configurations report only the parent job and
    # .batch/.extern steps. Use the parent final state in that case,
    # never the .extern state; no multi-task array uses this fallback.
    if "_" in str(jobid) and str(jobid).endswith("_1") and parent in entries:
        return ("terminal",*entries[parent])
    return ("unknown","NO_ACCOUNTING",None)


def is_remote_403_failure(jobid, root):
    """Read only the current job log, never classify prior runs' logs."""
    log=root/"logs"/("lr_stream_"+str(jobid)+".log")
    try:
        with log.open("rb") as fh:
            fh.seek(0,2)
            fh.seek(max(0,fh.tell()-131072))
            tail=fh.read().decode("utf-8",errors="replace")
    except OSError:
        return False
    return ("status=403" in tail.replace(" ", "") or "Forbidden" in tail)


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



def merge_complete(group, rows, out_root):
    sample,platform=group
    target=out_root/sample/platform/"numt_target_merged.bam"
    audit=out_root/sample/platform/"numt_target_merged.audit.json"
    if not (isfile(target) and isfile(Path(str(target)+".bai")) and isfile(audit)):
        return False
    a=json.loads(audit.read_text())
    wanted={r["run"] for r in rows}
    return (a.get("sample")==sample and a.get("platform")==platform
            and a.get("partial") is False
            and set(a.get("included_runs",[]))==wanted
            and a.get("n_included_runs")==len(wanted)
            and a.get("target_record_counts_verified") is True
            and set(a.get("per_run_target_records",{}))==wanted)


def submit_merge(group, root, manifest_path, state_dir, clean):
    sample,platform=group
    batch=state_dir/"jobs"/("merge_"+sample+"_"+platform+".tsv")
    with batch.open("w") as f:
        f.write("sample\tplatform\n"+sample+"\t"+platform+"\n")
    output=cmd([
        "sbatch","--parsable","--job-name=ln_mg_"+sample+"_"+platform,
        "--array=1-1",
        "--export=ALL,LR_PROJECT_ROOT="+str(root)+
        ",NUMT_MERGE_BATCH="+str(batch)+
        ",NUMT_MERGE_MANIFESTS="+str(manifest_path)+
        ",NUMT_MERGE_OUT_ROOT="+str(root/"results"/"sample_platform_bams")+
        ",NUMT_MERGE_ALLOW_PARTIAL=0,NUMT_MERGE_CLEANUP="+("1" if clean else "0")+
        ",NUMT_MERGE_MIN_CANDIDATES=10",
        str(root/"slurm"/"merge_sample_platform_numt.sbatch")
    ])
    parent=output.split(";",1)[0]
    if not parent.isdigit():
        raise RuntimeError("Bad sbatch job id: "+output)
    return parent+"_1"


def check_queue_group(group):
    sample,platform=group
    name="ln_mg_"+sample+"_"+platform
    data=cmd(["squeue","-h","-u",os.environ["USER"],"-o","%j|%i"])
    for line in data.splitlines():
        title,_,jid=line.partition("|")
        if title==name:
            return jid
    return None

def gib_free(path):
    s=os.statvfs(path)
    return s.f_bavail*s.f_frsize/(1024**3)


def run_once(a,do_submit,root,manifest,state_dir,manifest_path):
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
        manifest_hash=hashlib.sha256(json.dumps(manifest,sort_keys=True).encode()).hexdigest()
        previous=state.get("manifest_sha256")
        if previous and previous!=manifest_hash:
            raise RuntimeError("Manifest changed since controller start; stop controller and reconcile saved state before resuming")
        state["manifest_sha256"]=manifest_hash
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
                    observed_jobid=info["jobid"]
                    kind,slurm,exitcode=query_slurm(observed_jobid)
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
                        elif stage=="align" and exitcode==75 and is_remote_403_failure(observed_jobid,root):
                            info["status"]="BLOCKED_REMOTE_403"
                            info["last_error"]="ENA_HTTP_403: verify FTP URL / ENA access before retry"
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
        # Once every run in a biological sample/platform is fully validated,
        # archive its NUMT + chrM alignments with provenance. By design,
        # this is NOT triggered on partial run groups.
        merge_root=root/"results"/"sample_platform_bams"
        active_merges=0
        for group in order:
            key="merge:"+group[0]+":"+group[1]
            info=jobs.setdefault(key,{"attempts":0,"status":"PENDING"})
            finished=merge_complete(group,groups[group],merge_root)
            target=merge_root/group[0]/group[1]/"numt_target_merged.bam"
            if finished:
                info["status"]="DONE"
            elif target.exists() or (target.parent/"numt_target_merged.audit.json").exists():
                info["status"]="BLOCKED_EXISTING_INCOMPLETE_ARCHIVE"
            elif info.get("jobid"):
                kind,slurm,exitcode=query_slurm(info["jobid"])
                if kind=="active":
                    info["status"]="ACTIVE"
                    active_merges+=1
                elif kind=="unknown":
                    info["status"]="WAIT_ACCOUNTING"
                    active_merges+=1
                else:
                    info.pop("jobid")
                    if slurm=="COMPLETED":
                        info["status"]="BLOCKED_OUTPUT_MISSING"
                    elif slurm in RETRY_STATES:
                        info["status"]="RETRYABLE"
                        info["next_retry_time"]=time.time()+a.retry_delay_minutes*60
                    else:
                        info["status"]="BLOCKED_FAILURE"
                        info["last_error"]=slurm+":"+str(exitcode)
            elif info["status"]=="SUBMITTING":
                found=check_queue_group(group)
                if found:
                    info.update(jobid=found,status="ACTIVE",
                                attempts=info["attempts"]+1)
                    active_merges+=1
                else:
                    info["status"]="WAIT_SUBMISSION_RECONCILIATION"
            if info["status"] not in ("PENDING","RETRYABLE"):
                continue
            if not all(output_complete(r,"validate",results) for r in groups[group]):
                continue
            if active_merges>=a.max_merge or not disk_ok:
                continue
            if info["attempts"]>=a.max_attempts:
                info["status"]="MAX_ATTEMPTS"
                continue
            if info.get("next_retry_time",0)>time.time():
                continue
            old=check_queue_group(group)
            if old:
                info.update(jobid=old,status="ACTIVE",
                            external_discovered=True)
                active_merges+=1
                continue
            print(("SUBMIT" if do_submit else "WOULD_SUBMIT"),
                  "merge",group,flush=True)
            if do_submit:
                info["status"]="SUBMITTING"
                atomic_write(state_path,state)
                jid=submit_merge(group,root,manifest_path,state_dir,a.cleanup_if_validated)
                info.update(jobid=jid,status="ACTIVE",
                            attempts=info["attempts"]+1,submitted_utc=stamp())
                atomic_write(state_path,state)
            active_merges+=1

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
        all_done=all(jobs[status_key("validate",r["run"])]["status"]=="DONE" for r in tasks) and all(jobs["merge:"+g[0]+":"+g[1]]["status"]=="DONE" for g in order)
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
    ap.add_argument("--max-merge",type=int,default=1)
    ap.add_argument("--cleanup-if-validated",action="store_true",help="Allow irreversible full BAM deletion only after strict >90% gate")
    ap.add_argument("--max-attempts",type=int,default=5)
    ap.add_argument("--retry-delay-minutes",type=int,default=15)
    ap.add_argument("--submit",action="store_true",help="Allow Slurm submissions; default preview")
    ap.add_argument("--watch",action="store_true",help="Poll repeatedly; run via Slurm controller job")
    ap.add_argument("--poll-seconds",type=int,default=180)
    ap.add_argument("--watch-hours",type=float,default=18)
    a=ap.parse_args()
    if not all(x>0 for x in (a.min_free_gib,a.max_align,a.max_validate,
                              a.max_merge,a.max_attempts,a.retry_delay_minutes,
                              a.poll_seconds,a.watch_hours)):
        ap.error("All resource and timing parameters must be >0")
    if a.watch and not a.submit:
        ap.error("--watch requires --submit")
    if a.submit:
        old_names={"lr_numt_stream","primate_numt_validation"}
        active=cmd(["squeue","-h","-u",os.environ["USER"],"-o","%j"])
        if old_names & set(active.splitlines()):
            raise SystemExit("Existing legacy long-read Slurm jobs detected. Finish/cancel them before starting cohort manager to prevent duplicate submissions.")
    root=ROOT
    manifest_path=(root/a.manifest).resolve()
    manifest=load_manifest([manifest_path])
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
            finished=run_once(a,a.submit,root,manifest,state_dir,manifest_path)
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
