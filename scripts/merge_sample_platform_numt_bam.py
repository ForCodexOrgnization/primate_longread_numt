#!/usr/bin/env python3
"""Make one coordinate-sorted indexed NUMT-target BAM per biological sample/platform.

Reads are selected from independently aligned, indexed run BAMs using the
UNION of candidate nuclear NUMT intervals (+padding) and entire chrM.
No FASTQ is read; no per-run BAM is modified or deleted.

Every selected read receives @RG ID=<run>, SM=<sample>, PL=<platform>
and an RG tag, preserving origin for IGV and post-hoc grouping. This creates
a review/archive BAM, NOT a biological sample-level de novo NUMT call.
"""
import argparse
import csv
import fcntl
import json
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def run_cmd(cmd, dry=False):
    print("+", " ".join(map(str, cmd)), flush=True)
    if not dry:
        subprocess.run([str(x) for x in cmd], check=True)


def check_nonempty(path):
    p = Path(path)
    return p.is_file() and p.stat().st_size > 0


def merge_intervals(intervals):
    """Union only intervals on the same chromosome, never across contigs."""
    merged = []
    for chrom, start, end in sorted(intervals):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            merged[-1] = (chrom, merged[-1][1], max(end, merged[-1][2]))
        else:
            merged.append((chrom, start, end))
    return merged


def gather(rows, root, padding, chrm, allow_partial):
    """Validate inputs before writing anything. Returns selected and missing."""
    available, missing, nuclear = [], [], []
    references = {r["reference"] for r in rows}
    if len(references) != 1:
        raise ValueError("Reference path mismatch across runs; refusing BAM merge")
    besthits = {r["besthit"] for r in rows}
    if len(besthits) != 1:
        raise ValueError("Besthit candidate source differs across runs; inspect before merging")
    for r in rows:
        base = root / r["sample"] / r["run"]
        bam = base / "aligned.bam"
        bai = base / "aligned.bam.bai"
        loci = base / "chrm_candidate_validation" / "chrM_candidate_locus_validation.tsv"
        required = [bam, bai, base / "complete.flag", loci]
        absent = [str(p) for p in required if not (p.is_file() if p.name == "complete.flag" else check_nonempty(p))]
        if absent:
            missing.append({"run":r["run"], "reason":";".join(absent)})
            continue
        candidates = read_tsv(loci)
        if not candidates:
            missing.append({"run":r["run"], "reason":"empty_locus_table"})
            continue
        for loc in candidates:
            nuclear.append((loc["nuclear_chrom"], int(loc["nuclear_start0"]), int(loc["nuclear_end0"])))
        available.append((r, bam))
    if missing and not allow_partial:
        raise RuntimeError("Some manifest runs are not ready: " + json.dumps(missing[:15]) +
                           ". Use --allow-partial only for exploratory output.")
    if not available:
        raise RuntimeError("No aligned runs with NUMT candidate evidence available")
    return available, missing, nuclear


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", action="append", required=True, help="Repeat for pilot and expanded manifests")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--platform", choices=("ONT", "PACBIO"), required=True)
    ap.add_argument("--results-root", default="results")
    ap.add_argument("--out-root", default="results/sample_platform_bams")
    ap.add_argument("--padding", type=int, default=1000)
    ap.add_argument("--chrM-contig", default="chrM")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--allow-partial", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    if args.padding < 100 or args.threads < 1:
        ap.error("padding >=100 and threads >=1 required")
    if not re.fullmatch(r"[A-Za-z0-9._-]+", args.sample):
        ap.error("Invalid sample identifier")
    if not shutil.which("samtools"):
        ap.error("samtools not found in PATH")
    for manifest in args.manifest:
        if not Path(manifest).is_file():
            ap.error(f"Missing manifest: {manifest}")
    selected = {}
    for path in args.manifest:
        for row in read_tsv(path):
            if row["sample"] != args.sample or row["platform"] != args.platform:
                continue
            run = row["run"]
            if not re.fullmatch(r"[A-Za-z0-9._-]+", run):
                raise ValueError(f"Invalid run accession: {run}")
            if run in selected and selected[run] != row:
                # Same run from two manifests may differ only in FASTQ URL order;
                # still do not allow ambiguous entries.
                raise ValueError(f"Conflicting manifest entries for {run}")
            selected[run] = row
    if not selected:
        ap.error("No entries found for sample/platform")
    rows = list(selected.values())
    species = {x["species"] for x in rows}
    if len(species) != 1:
        raise ValueError("Different species across runs")
    root = Path(args.results_root).resolve()
    ready, missing, nuclear = gather(rows, root, args.padding, args.chrM_contig, args.allow_partial)
    merged_loci = merge_intervals((c,max(0,s-args.padding),e+args.padding)
                                  for c,s,e in nuclear)
    destdir = Path(args.out_root).resolve() / args.sample / args.platform
    dest = destdir / "numt_target_merged.bam"
    audit = destdir / "numt_target_merged.audit.json"
    report = dict(sample=args.sample, species=next(iter(species)), platform=args.platform,
                  n_manifest_runs=len(rows), n_included_runs=len(ready),
                  partial=bool(missing), missing_runs=missing,
                  included_runs=[r["run"] for r,_ in ready],
                  reference=rows[0]["reference"], besthit=rows[0]["besthit"],
                  padding=args.padding, chrM_contig=args.chrM_contig,
                  n_merged_nuclear_intervals=len(merged_loci),
                  source_bams=[str(b) for _,b in ready],
                  original_run_bams_retained=True, fastq_files_untouched=True)
    print(json.dumps(report, indent=2))
    if args.dry_run:
        print("DRY RUN: no BAM generated")
        return
    destdir.mkdir(parents=True,exist_ok=True)
    lock = destdir / ".merge.lock"
    with lock.open("w") as guard:
        fcntl.flock(guard, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if dest.exists() or audit.exists():
            raise RuntimeError(f"Output exists; refusing overwrite: {dest}. Archive/remove old output explicitly.")
        with tempfile.TemporaryDirectory(prefix=".merge_work_",dir=destdir) as td:
            work = Path(td)
            chrom_lengths = None
            for r,bam in ready:
                # Compare @SQ headers, including contig names, lengths, order.
                raw = subprocess.check_output(["samtools","view","-H",str(bam)],text=True)
                sq = [line for line in raw.splitlines() if line.startswith("@SQ\t")]
                if chrom_lengths is None:
                    chrom_lengths = sq
                elif sq != chrom_lengths:
                    raise ValueError(f"Reference header / contig mismatch for {r['run']}")
            sq_map = {}
            for line in chrom_lengths:
                attrs = dict(x.split(":",1) for x in line.split("\t")[1:])
                sq_map[attrs["SN"]] = int(attrs["LN"])
            if args.chrM_contig not in sq_map:
                raise ValueError("chrM contig missing from BAM header")
            intervals = merge_intervals(
                [(chrom,start,min(end,sq_map.get(chrom,0))) for chrom,start,end in merged_loci]
                + [(args.chrM_contig,0,sq_map[args.chrM_contig])])
            if any(chrom not in sq_map or end <= start for chrom,start,end in intervals):
                raise ValueError("Candidate interval not present in reference or invalid coordinates")
            bed=work/"regions.bed"
            with bed.open("w") as fh:
                for chrom,start,end in intervals:
                    fh.write(f"{chrom}\t{start}\t{end}\n")
            inputs=[]
            for i,(r,bam) in enumerate(ready):
                extracted=work/f"{i:05d}.extracted.bam"
                tagged=work/f"{i:05d}.tagged.bam"
                run_cmd(["samtools","quickcheck",str(bam)])
                run_cmd(["samtools","view","-b","-L",str(bed),"-o",str(extracted),str(bam)])
                run_cmd(["samtools","addreplacerg","-m","overwrite_all",
                         "-r",f"ID:{r['run']}","-r",f"SM:{args.sample}",
                         "-r",f"PL:{args.platform}",
                         "-o",str(tagged),str(extracted)])
                extracted.unlink()
                inputs.append(tagged)
            inlist=work/"bams.list"
            inlist.write_text("".join(str(x)+"\n" for x in inputs))
            merged=work/"merged.bam"
            # All per-run files are coordinate sorted and share reference SQ;
            # samtools merge keeps coordinate ordering and distinct RG IDs.
            run_cmd(["samtools","merge","-@",
                     str(args.threads),"-b",str(inlist),"-o",str(merged)])
            run_cmd(["samtools","quickcheck",str(merged)])
            hdr = subprocess.check_output(["samtools","view","-H",str(merged)],text=True)
            rg_ids = set()
            for line in hdr.splitlines():
                if line.startswith("@RG\\t"):
                    fields = dict(x.split(":",1) for x in line.split("\\t")[1:])
                    rg_ids.add(fields.get("ID"))
            expected_rg = {r["run"] for r,_ in ready}
            if not expected_rg.issubset(rg_ids):
                raise RuntimeError(f"Merged BAM missing run read groups: {sorted(expected_rg-rg_ids)}")
            run_cmd(["samtools","index","-@",str(args.threads),str(merged)])
            run_cmd(["samtools","idxstats",str(merged)])
            if not check_nonempty(merged) or not check_nonempty(Path(str(merged)+".bai")):
                raise RuntimeError("Merged BAM/index verification failed")
            report["n_target_intervals"]=len(intervals)
            report["output_bytes"]=merged.stat().st_size
            # Promote data only after verified merge/index; preserve original BAMs.
            os.replace(merged, dest)
            os.replace(Path(str(merged)+".bai"),Path(str(dest)+".bai"))
            audit.write_text(json.dumps(report,indent=2)+"\n")
        print("COMPLETE:",dest)
        print("Saved:",audit)
        print("Per-run full-genome BAMs were NOT deleted.")


if __name__ == "__main__":
    main()
