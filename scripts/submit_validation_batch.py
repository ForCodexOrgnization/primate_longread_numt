#!/usr/bin/env python3
"""Submit only fully aligned pilot runs for two-stage NUMT validation."""
import argparse
import csv
import datetime
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def nonempty(file):
    file = Path(file)
    return file.is_file() and file.stat().st_size > 0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", default="data/pilot_remote_manifest.tsv")
    p.add_argument("--batch-size", type=int, default=5)
    p.add_argument("--concurrent", type=int, default=2)
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()
    if min(a.batch_size, a.concurrent) < 1:
        p.error("batch-size and concurrent must be positive")
    source = Path(a.manifest)
    if not source.is_absolute():
        source = ROOT / source
    with source.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    queued, waiting, done = [], [], []
    for r in rows:
        base = ROOT / "results" / r["sample"] / r["run"]
        required = [base / "complete.flag", base / "aligned.bam",
                    base / "aligned.bam.bai", r["reference"],
                    str(r["reference"]) + ".fai", r["besthit"]]
        if not all(nonempty(f) or (Path(f).name == "complete.flag" and Path(f).is_file()) for f in required):
            waiting.append((r["sample"], r["run"]))
            continue
        n = base / "chrm_candidate_validation"
        c = base / "longread_chrM_homology"
        outputs = [n / "chrM_candidate_locus_validation.tsv",
                   n / "chrM_candidate_read_evidence.tsv",
                   c / "per_locus_chrM_homology.tsv",
                   c / "per_read_chrM_homology.tsv"]
        if all(nonempty(f) for f in outputs):
            done.append((r["sample"], r["run"]))
        elif len(queued) < a.batch_size:
            queued.append(r)
    print(f"Manifest runs={len(rows)}, aligned-not-ready={len(waiting)}, validated={len(done)}, submit={len(queued)}")
    if not queued:
        return
    batches = ROOT / "data" / "batches"
    batches.mkdir(parents=True, exist_ok=True)
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    dest = batches / f"validation_{stamp}.tsv"
    with dest.open("w", newline="") as handle:
        wr = csv.DictWriter(handle, fieldnames=list(queued[0]), delimiter="\t", lineterminator="\n")
        wr.writeheader()
        wr.writerows(queued)
    command = ["sbatch", f"--array=1-{len(queued)}%{min(a.concurrent, len(queued))}",
               f"--export=ALL,LR_PROJECT_ROOT={ROOT},VALIDATION_MANIFEST={dest}",
               str(ROOT / "slurm" / "validate_pilot_array.sbatch")]
    print(" ".join(command))
    if not a.dry_run:
        subprocess.run(command, check=True, cwd=ROOT)
    else:
        print("DRY RUN: no jobs submitted")


if __name__ == "__main__":
    main()
