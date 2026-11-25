#!/usr/bin/env python3
"""Audit all ENA-resolved cohort runs, matched Illumina, and long-read groups.

Read-only; produces reports/planning lists, never queues jobs. Uses all resolved
cohort runs (including Illumina and pilot), not the non-pilot-only file.
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path


def load(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fastq_total(value):
    """ENA fastq_bytes is semicolon-delimited per FASTQ; unknown != zero."""
    vals = [x.strip() for x in (value or "").split(";") if x.strip()]
    if not vals:
        return None
    try:
        sizes = [int(v) for v in vals]
    except ValueError:
        return None
    return sum(sizes) if all(v >= 0 for v in sizes) else None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--resolved", default="data/all_longread_resolved_runs.tsv")
    p.add_argument("--metadata", default="data/all_longread_samples.tsv")
    p.add_argument("--pilot-selection", default="data/pilot_run_selection.tsv")
    p.add_argument("--outdir", default="data/all_longread_audit")
    args = p.parse_args()
    out = Path(args.outdir)
    runs = load(args.resolved)
    samples = load(args.metadata)
    pilot = {r["run_accession"] for r in load(args.pilot_selection)}
    if not runs:
        raise SystemExit("ERROR: no resolved runs")
    samplemeta = {}
    for r in samples:
        k = r["input_accession"]
        if k in samplemeta:
            raise ValueError(f"Duplicate metadata sample: {k}")
        samplemeta[k] = r
    bysample = defaultdict(list)
    bygroup = defaultdict(list)
    problems = []
    seen = {}
    for r in runs:
        accession = r["run_accession"]
        sid = r["sample"]
        key = (sid, r["species"])
        if accession in seen:
            problems.append(dict(run=accession, sample=sid, reason="duplicate_run", detail=str(seen[accession])))
            continue
        seen[accession] = key
        if sid not in samplemeta:
            problems.append(dict(run=accession, sample=sid, reason="missing_metadata", detail=""))
        elif r["species"] != samplemeta[sid]["pipeline_species"]:
            problems.append(dict(run=accession, sample=sid, reason="species_mismatch", detail=r["species"]))
        platform = r.get("run_platform","").upper()
        kind = "ONT" if platform == "OXFORD_NANOPORE" else "PACBIO" if "PACBIO" in platform else "ILLUMINA" if platform == "ILLUMINA" else "OTHER"
        r["_kind"] = kind
        r["_bytes"] = fastq_total(r.get("fastq_bytes",""))
        if kind in ("ONT","PACBIO") and r.get("read_type") != kind:
            problems.append(dict(run=accession,sample=sid,reason="read_type_mismatch",detail=r.get("read_type","")))
        if r.get("lookup_error"):
            problems.append(dict(run=accession,sample=sid,reason="ena_lookup_error",detail=r["lookup_error"]))
        bysample[sid].append(r)
        if kind in ("ONT", "PACBIO"):
            bygroup[(sid, r["species"], kind)].append(r)

    sample_rows=[]
    for sid, meta in sorted(samplemeta.items()):
        rr=bysample.get(sid,[])
        illum=[r for r in rr if r["_kind"]=="ILLUMINA"]
        illumfq=[r for r in illum if r.get("fastq_ftp","").strip()]
        longread=[r for r in rr if r["_kind"] in ("ONT","PACBIO")]
        mplatform=meta.get("sequencing_platform","").upper()
        claimed="ILLUMINA" in mplatform
        status=("VERIFIED_ILLUMINA_FASTQ" if illumfq else
                "ILLUMINA_RUN_FASTQ_MISSING" if illum else
                "METADATA_ONLY_ILLUMINA" if claimed else "NO_ILLUMINA_EVIDENCE")
        sample_rows.append(dict(sample=sid,biosample=meta.get("biosample_accession",""),
                                species=meta["pipeline_species"],metadata_illumina=claimed,
                                short_read_status=status,n_illumina_runs=len(illum),
                                n_illumina_fastq_runs=len(illumfq),
                                illumina_runs=";".join(r["run_accession"] for r in illumfq),
                                n_ont_runs=sum(r["_kind"]=="ONT" for r in longread),
                                n_pacbio_runs=sum(r["_kind"]=="PACBIO" for r in longread),
                                n_longread_fastq_runs=sum(bool(r.get("fastq_ftp","").strip()) for r in longread)))
    save(out/"sample_shortread_status.tsv",sample_rows,
         ["sample","biosample","species","metadata_illumina","short_read_status",
          "n_illumina_runs","n_illumina_fastq_runs","illumina_runs",
          "n_ont_runs","n_pacbio_runs","n_longread_fastq_runs"])

    groups=[]
    planning=[]
    sr={x["sample"]:x["short_read_status"] for x in sample_rows}
    for (sid,species,kind), rr in sorted(bygroup.items()):
        available=[r for r in rr if r.get("fastq_ftp","").strip() and not r.get("lookup_error")]
        nonpilot=[r for r in available if r["run_accession"] not in pilot]
        sized=[r for r in available if r["_bytes"] is not None]
        total_bytes=sum(r["_bytes"] for r in sized)
        groups.append(dict(sample=sid,species=species,platform=kind,
                           short_read_status=sr.get(sid,"MISSING_METADATA"),
                           n_runs=len(rr),n_fastq_runs=len(available),
                           n_pilot_runs=sum(r["run_accession"] in pilot for r in rr),
                           n_nonpilot_fastq_runs=len(nonpilot),
                           n_size_known=len(sized),
                           total_fastq_gb_known=round(total_bytes/1e9,3),
                           missing_size_runs=len(available)-len(sized),
                           run_ids=";".join(r["run_accession"] for r in rr)))
        # Size is a download/storage proxy, NOT read quality or coverage.
        ordered=sorted(nonpilot,key=lambda r: (r["_bytes"] is None, -(r["_bytes"] or 0),r["run_accession"]))
        for index,r in enumerate(ordered,1):
            planning.append(dict(sample=sid,species=species,platform=kind,
                                 run=r["run_accession"],rank_by_fastq_size=index,
                                 fastq_gb=round(r["_bytes"]/1e9,3) if r["_bytes"] is not None else "",
                                 size_known=r["_bytes"] is not None,
                                 short_read_status=sr.get(sid,"MISSING_METADATA"),
                                 run_model=r.get("run_model",""),preset_review_required=(kind=="PACBIO"),
                                 note="Size rank only; PacBio CLR/HiFi must be verified"))
    save(out/"sample_platform_groups.tsv",groups,
         ["sample","species","platform","short_read_status","n_runs","n_fastq_runs",
          "n_pilot_runs","n_nonpilot_fastq_runs","n_size_known","total_fastq_gb_known",
          "missing_size_runs","run_ids"])
    save(out/"nonpilot_run_priority_review.tsv",planning,
         ["sample","species","platform","run","rank_by_fastq_size","fastq_gb",
          "size_known","short_read_status","run_model","preset_review_required","note"])
    save(out/"audit_problems.tsv",problems,["run","sample","reason","detail"])
    from collections import Counter
    print("Metadata biological samples:",len(sample_rows))
    print("Resolved unique ENA runs:",len(seen))
    print("Short-read status:",dict(Counter(x["short_read_status"] for x in sample_rows)))
    print("Sample-platform groups:",len(groups))
    print("Non-pilot FASTQ runs to review:",len(planning))
    print("Problem flags:",dict(Counter(x["reason"] for x in problems)))
    print("Reports:",out)
    print("No jobs submitted. A FASTQ size ranking is NOT evidence of high coverage or quality.")


if __name__ == "__main__":
    main()
