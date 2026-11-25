#!/usr/bin/env python3
"""Run post-alignment nuclear and chrM NUMT checks for one pilot manifest row."""
import argparse
import csv
import subprocess
import sys
from pathlib import Path


def nonempty(path):
    path = Path(path)
    return path.is_file() and path.stat().st_size > 0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True)
    p.add_argument("--task-id", type=int, required=True, help="1-based data-row index")
    p.add_argument("--root", default=str(Path(__file__).resolve().parent.parent))
    args = p.parse_args()
    root = Path(args.root).resolve()
    with open(args.manifest, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not 1 <= args.task_id <= len(rows):
        p.error(f"task-id {args.task_id} outside manifest (1..{len(rows)})")
    row = rows[args.task_id - 1]
    sample, run = row["sample"], row["run"]
    base = root / "results" / sample / run
    bam = base / "aligned.bam"
    ref = Path(row["reference"])
    besthit = Path(row["besthit"])
    for name in (bam, Path(str(bam) + ".bai"), ref, Path(str(ref) + ".fai"), besthit, base / "complete.flag"):
        if not nonempty(name) and not (name.name == "complete.flag" and name.is_file()):
            raise SystemExit(f"[PRECHECK] Missing/incomplete input: {name}")
    nuclear = base / "chrm_candidate_validation"
    homology = base / "longread_chrM_homology"
    nuclear_files = (nuclear / "chrM_candidate_locus_validation.tsv",
                     nuclear / "chrM_candidate_read_evidence.tsv")
    chrM_files = (homology / "per_locus_chrM_homology.tsv",
                  homology / "per_read_chrM_homology.tsv")
    if not all(nonempty(f) for f in nuclear_files):
        cmd = [sys.executable, str(root / "scripts/validate_chrm_candidates.py"),
               "--bam", str(bam), "--reference", str(ref),
               "--besthit", str(besthit), "--outdir", str(nuclear)]
        print(f"[NUCLEAR] {sample} {run} {row['species']} {row['preset']}", flush=True)
        subprocess.run(cmd, check=True)
    else:
        print(f"[SKIP NUCLEAR] {sample} {run}", flush=True)
    if not all(nonempty(f) for f in chrM_files):
        cmd = [sys.executable, str(root / "scripts/validate_read_chrm_homology.py"),
               "--bam", str(bam), "--reference", str(ref),
               "--nuclear-validation", str(nuclear), "--outdir", str(homology),
               "--threads", "2"]
        print(f"[CHRM] {sample} {run}", flush=True)
        subprocess.run(cmd, check=True)
    else:
        print(f"[SKIP CHRM] {sample} {run}", flush=True)
    print(f"[PILOT VALIDATION COMPLETE] {sample} {run}", flush=True)


if __name__ == "__main__":
    main()
