#!/usr/bin/env python3
"""Sample-level NUMT validation by pooling distinct read evidence across runs.

Input: one or more eight-column remote manifests. Uses completed per-run
nuclear and chrM evidence, never changes the run-level outputs or BAMs.

Strict classification counts unique (run, read) pairs passing ALL nuclear
criteria, and the SUBSET of those same pairs passing expected chrM homology.
Reads from different ONT and PacBio platforms are also summarized separately.

Default thresholds exactly match validate_chrm_candidates.py:
MAPQ >=20; identity >=.80; aligned fraction >=.90; both 100-bp
flanks >=80 aligned bases; max internal deletion <=20; >=2 reads.
Only originally retained per-read evidence can be re-thresholded. Reads with
MAPQ<20 and chrM BLAST results for non-original-strict reads are unavailable.
"""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def load(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".partial")
    with tmp.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


def num(row, key, default=-1):
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return default


def signature(row):
    # Locus IDs can be reused or change across different short-read samples.
    # Require exact reference coordinates and chrM target for safe merging.
    return (row["nuclear_chrom"], row["nuclear_start0"],
            row["nuclear_end0"], row["chrM_hit"], row["chrM_start"],
            row["chrM_end"], row["group"])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", action="append", required=True,
                    help="May repeat for pilot + nonpilot manifest")
    ap.add_argument("--results-root", default="results")
    ap.add_argument("--outdir", default="results/sample_level_validation")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-identity", type=float, default=.80)
    ap.add_argument("--min-locus-fraction", type=float, default=.90)
    ap.add_argument("--min-flank-aligned", type=int, default=80)
    ap.add_argument("--max-deletion", type=int, default=20)
    ap.add_argument("--min-support", type=int, default=2)
    args = ap.parse_args()
    if args.min_mapq < 20:
        ap.error("MAPQ<20 cannot be recovered from stored per-read evidence")
    if not (0 < args.min_identity <= 1 and 0 < args.min_locus_fraction <= 1):
        ap.error("fraction thresholds must be within (0,1]")
    if not (0 <= args.min_flank_aligned <= 100 and args.max_deletion >= 0 and args.min_support >= 1):
        ap.error("invalid flank/deletion/support settings")

    root, out = Path(args.results_root), Path(args.outdir)
    manifest = {}
    for file in args.manifest:
        for r in load(file):
            key = (r["sample"], r["run"])
            if key in manifest:
                old = manifest[key]
                for field in ("species", "platform", "reference", "besthit"):
                    if old[field] != r[field]:
                        raise ValueError(f"Conflicting manifest for {key}, field={field}")
            else:
                manifest[key] = r
    if not manifest:
        raise SystemExit("No manifest entries")

    records = defaultdict(lambda: dict(runs=set(), platforms=set(), nuclear=set(),
                                       chrM=set(), nuclear_by_platform=defaultdict(set),
                                       chrm_by_platform=defaultdict(set)))
    sample_runs = defaultdict(set)
    sample_platforms = defaultdict(set)
    status = []
    reference_by_sample = {}
    # Read-level evidence key: (run, read name). This prevents collisions
    # between independent run accessions, and deduplicates repeated records.
    for (sample, run), mr in sorted(manifest.items()):
        species, platform = mr["species"], mr["platform"]
        s = (sample, species)
        sample_runs[s].add(run)
        sample_platforms[s].add(platform)
        ref = mr["reference"]
        if s in reference_by_sample and reference_by_sample[s] != ref:
            raise ValueError(f"Reference mismatch across runs for {sample}: {ref}")
        reference_by_sample[s] = ref
        base = root / sample / run
        nbase = base / "chrm_candidate_validation"
        cbase = base / "longread_chrM_homology"
        required = [nbase / "chrM_candidate_locus_validation.tsv",
                    nbase / "chrM_candidate_read_evidence.tsv",
                    cbase / "per_locus_chrM_homology.tsv",
                    cbase / "per_read_chrM_homology.tsv"]
        missing = [str(x) for x in required if not x.is_file() or x.stat().st_size == 0]
        if missing:
            status.append(dict(sample=sample,species=species,platform=platform,run=run,
                               status="PENDING",detail=";".join(missing)))
            continue
        locs = load(required[0])
        reads = load(required[1])
        clocs = load(required[2])
        chrreads = load(required[3])
        li = {x["locus_id"]: x for x in locs}
        ci = {x["locus_id"]: x for x in clocs}
        if len(li) != len(locs) or len(ci) != len(clocs) or set(li) != set(ci):
            raise ValueError(f"Inconsistent locus IDs in {sample}/{run}")
        # This legacy output tests at most 40 strict nuclear reads per locus.
        chr_eligible = defaultdict(set)
        for x in chrreads:
            if x["locus_id"] not in li:
                raise ValueError(f"Unexpected chrM read locus in {sample}/{run}")
            if str(x.get("expected_chrM_homology", "0")) == "1":
                chr_eligible[x["locus_id"]].add((run, x["read"]))
        nuclear_eligible = defaultdict(set)
        for x in reads:
            if x["locus_id"] not in li:
                raise ValueError(f"Unexpected nuclear read locus in {sample}/{run}")
            if (num(x,"mapq") >= args.min_mapq
                and num(x,"locus_identity_to_nuclear_ref") >= args.min_identity
                and num(x,"locus_fraction") >= args.min_locus_fraction
                and num(x,"left_aligned") >= args.min_flank_aligned
                and num(x,"right_aligned") >= args.min_flank_aligned
                and num(x,"max_internal_deletion", 1e9) <= args.max_deletion):
                nuclear_eligible[x["locus_id"]].add((run, x["read"]))
        for lid, lr in li.items():
            key = (sample, species, lid)
            sig = signature(lr)
            if key in records and records[key].get("signature") != sig:
                raise ValueError(f"Locus discrepancy between runs: {key}")
            record = records[key]
            record["signature"] = sig
            record["runs"].add(run)
            record["platforms"].add(platform)
            nuclei = nuclear_eligible.get(lid, set())
            # chrM read must be independently nuclear-eligible under selected
            # nuclear settings. An expected chrM match alone is NOT sufficient.
            both = chr_eligible.get(lid, set()) & nuclei
            record["nuclear"].update(nuclei)
            record["chrM"].update(both)
            record["nuclear_by_platform"][platform].update(nuclei)
            record["chrm_by_platform"][platform].update(both)
        status.append(dict(sample=sample,species=species,platform=platform,run=run,
                           status="VALIDATED",detail=f"{len(li)} loci"))

    if not records:
        raise SystemExit("No completed per-run evidence in requested manifests")
    locus_out=[]
    for (sample,species,lid),x in sorted(records.items()):
        n_nuclear=len(x["nuclear"])
        n_both=len(x["chrM"])
        if n_both > n_nuclear:
            raise AssertionError("chrM evidence exceeds nuclear evidence")
        row=dict(sample=sample,species=species,locus_id=lid,
                 nuclear_chrom=x["signature"][0],nuclear_start0=x["signature"][1],
                 nuclear_end0=x["signature"][2],chrM_hit=x["signature"][3],
                 chrM_start=x["signature"][4],chrM_end=x["signature"][5],
                 group=x["signature"][6],n_validated_runs=len(x["runs"]),
                 platforms=";".join(sorted(x["platforms"])),
                 n_nuclear_strict_reads=n_nuclear,
                 n_expected_chrM_reads=n_both,
                 nuclear_status="SUPPORTED" if n_nuclear>=args.min_support else "NOT_SUPPORTED",
                 chrM_status="EXPECTED_REGION_SUPPORTED" if n_both>=args.min_support else "NOT_SUPPORTED",
                 final_status="VALIDATED" if n_both>=args.min_support else "NOT_VALIDATED")
        locus_out.append(row)

    locus_fields=["sample","species","locus_id","nuclear_chrom","nuclear_start0",
                  "nuclear_end0","chrM_hit","chrM_start","chrM_end","group",
                  "n_validated_runs","platforms","n_nuclear_strict_reads",
                  "n_expected_chrM_reads","nuclear_status","chrM_status","final_status"]
    write(out/"sample_locus_validation.tsv",locus_out,locus_fields)
    platform_out=[]
    for (sample,species,lid),x in sorted(records.items()):
        for platform in sorted(x["platforms"]):
            n=len(x["nuclear_by_platform"][platform])
            m=len(x["chrm_by_platform"][platform])
            platform_out.append(dict(sample=sample,species=species,platform=platform,locus_id=lid,
                                     group=x["signature"][6],n_nuclear_strict_reads=n,
                                     n_expected_chrM_reads=m,
                                     final_status="VALIDATED" if m>=args.min_support else "NOT_VALIDATED"))
    write(out/"sample_platform_locus_validation.tsv",platform_out,
          ["sample","species","platform","locus_id","group","n_nuclear_strict_reads",
           "n_expected_chrM_reads","final_status"])
    summary=[]
    by_sample=defaultdict(list)
    for r in locus_out:
        by_sample[(r["sample"],r["species"])].append(r)
    for (sample,species),rs in sorted(by_sample.items()):
        summary.append(dict(sample=sample,species=species,
                            n_manifest_runs=len(sample_runs[(sample,species)]),
                            n_completed_runs=sum(z["sample"]==sample and z["species"]==species and z["status"]=="VALIDATED" for z in status),
                            platforms=";".join(sorted(sample_platforms[(sample,species)])),
                            n_loci=len(rs),
                            n_blast_qc_pass=sum(r["group"]=="BLAST_QC_PASS" for r in rs),
                            n_nuclear_supported=sum(r["nuclear_status"]=="SUPPORTED" for r in rs),
                            n_both_supported=sum(r["final_status"]=="VALIDATED" for r in rs)))
    write(out/"sample_summary.tsv",summary,
          ["sample","species","n_manifest_runs","n_completed_runs","platforms","n_loci",
           "n_blast_qc_pass","n_nuclear_supported","n_both_supported"])
    write(out/"run_status.tsv",status,["sample","species","platform","run","status","detail"])
    parameters=dict(min_mapq=args.min_mapq,min_identity=args.min_identity,
                    min_locus_fraction=args.min_locus_fraction,min_flank_aligned=args.min_flank_aligned,
                    max_deletion=args.max_deletion,min_support=args.min_support,
                    manifest=args.manifest,results_root=str(root),
                    limitation="chrM read homology was computed only for reads originally strict under run-level QC; max 40 per locus per run")
    (out/"parameters.json").write_text(json.dumps(parameters,indent=2)+"\n")
    print(f"Sample-level: {len(summary)} samples, {len(locus_out)} candidate loci; completed runs: {sum(x['status']=='VALIDATED' for x in status)}/{len(status)}")
    print(f"Validated loci: {sum(x['final_status']=='VALIDATED' for x in locus_out)}")
    print("Output:",out)
    print("NOTE: chrM homology for nuclear reads failing original strict criteria was not BLAST-tested. Rerun per-read BLAST before interpreting relaxed nuclear thresholds.")


if __name__=="__main__":
    main()
