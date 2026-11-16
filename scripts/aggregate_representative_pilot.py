#!/usr/bin/env python3
"""Aggregate representative-pilot status and per-run NUMT validation by species.

Each biological sample counts once regardless of the number of ONT/PacBio runs.
Technology comparisons stay separate and are not treated as independent samples.
"""
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def read(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def has(path):
    path = Path(path)
    return path.is_file() and path.stat().st_size > 0


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", default="data/pilot_remote_manifest.tsv")
    p.add_argument("--results-root", default="results")
    p.add_argument("--outdir", default="results/representative_pilot_summary")
    args = p.parse_args()
    root = Path(args.results_root)
    outdir = Path(args.outdir)
    run_rows = []
    locus_rows = []
    sample_runs = defaultdict(list)
    for r in read(args.manifest):
        sample, run, species = r["sample"], r["run"], r["species"]
        platform, preset = r["platform"], r["preset"]
        base = root / sample / run
        nuclear = base / "chrm_candidate_validation" / "chrM_candidate_locus_validation.tsv"
        homology = base / "longread_chrM_homology" / "per_locus_chrM_homology.tsv"
        aligned = has(base / "aligned.bam") and has(base / "aligned.bam.bai") and (base / "complete.flag").is_file()
        completed = aligned and has(nuclear) and has(homology)
        status = "VALIDATED" if completed else "ALIGNED_PENDING_VALIDATION" if aligned else "NOT_ALIGNED"
        sample_runs[(sample, species)].append((run, status, platform, preset))
        run_record = dict(sample=sample, species=species, run=run, platform=platform, preset=preset,
                          status=status, n_chrM_hit="", n_blast_qc_pass="", n_nuclear_supported="",
                          n_expected_chrM_supported="", n_both_supported="")
        if completed:
            nrows = read(nuclear)
            hrows = read(homology)
            ni = {x["locus_id"]: x for x in nrows}
            hi = {x["locus_id"]: x for x in hrows}
            if set(ni) != set(hi) or len(ni) != len(nrows) or len(hi) != len(hrows):
                raise ValueError(f"Locus mismatch for {sample}/{run}: nuclear={len(nrows)} chrM={len(hrows)}")
            run_record.update(
                n_chrM_hit=len(ni),
                n_blast_qc_pass=sum(x["group"] == "BLAST_QC_PASS" for x in nrows),
                n_nuclear_supported=sum(x["strict_status"] == "SUPPORTED" for x in nrows),
                n_expected_chrM_supported=sum(x["chrM_status"] == "EXPECTED_REGION_SUPPORTED" for x in hrows),
                n_both_supported=sum(ni[k]["strict_status"] == "SUPPORTED"
                                     and hi[k]["chrM_status"] == "EXPECTED_REGION_SUPPORTED" for k in ni))
            for key in sorted(ni):
                n, h = ni[key], hi[key]
                locus_rows.append(dict(sample=sample, species=species, run=run, platform=platform,
                                       preset=preset, locus_id=key, nuclear_chrom=n["nuclear_chrom"],
                                       nuclear_start0=n["nuclear_start0"], nuclear_end0=n["nuclear_end0"],
                                       group=n["group"], nuclear_strict_status=n["strict_status"],
                                       n_strict_sequence_support=n["n_strict_sequence_support"],
                                       chrM_status=h["chrM_status"],
                                       n_expected_chrM_homology=h["n_expected_chrM_homology"],
                                       both_supported=int(n["strict_status"] == "SUPPORTED"
                                                          and h["chrM_status"] == "EXPECTED_REGION_SUPPORTED")))
        run_rows.append(run_record)
    fields = ["sample", "species", "run", "platform", "preset", "status", "n_chrM_hit",
              "n_blast_qc_pass", "n_nuclear_supported", "n_expected_chrM_supported", "n_both_supported"]
    write(outdir / "run_summary.tsv", run_rows, fields)
    locus_fields = ["sample", "species", "run", "platform", "preset", "locus_id",
                    "nuclear_chrom", "nuclear_start0", "nuclear_end0", "group",
                    "nuclear_strict_status", "n_strict_sequence_support", "chrM_status",
                    "n_expected_chrM_homology", "both_supported"]
    write(outdir / "run_locus_evidence.tsv", locus_rows, locus_fields)
    loci_by_sample = defaultdict(lambda: defaultdict(list))
    for x in locus_rows:
        site = (x["nuclear_chrom"], x["nuclear_start0"], x["nuclear_end0"], x["group"])
        loci_by_sample[(x["sample"], x["species"])][site].append(x)
    samples = []
    for (sample, species), sruns in sorted(sample_runs.items()):
        loci = loci_by_sample[(sample, species)]
        complete_runs = [r for r, status, _, _ in sruns if status == "VALIDATED"]
        samples.append(dict(
            sample=sample, species=species, n_manifest_runs=len(sruns),
            n_validated_runs=len(complete_runs),
            platforms=";".join(sorted(set(p for _, _, p, _ in sruns))),
            validated_platforms=";".join(sorted(set(p for _, status, p, _ in sruns if status == "VALIDATED"))),
            n_loci_with_results=len(loci),
            n_any_nuclear_supported=sum(any(z["nuclear_strict_status"] == "SUPPORTED" for z in entries)
                                        for entries in loci.values()),
            n_any_chrM_supported=sum(any(z["chrM_status"] == "EXPECTED_REGION_SUPPORTED" for z in entries)
                                     for entries in loci.values()),
            n_any_both_supported=sum(any(int(z["both_supported"]) == 1 for z in entries)
                                     for entries in loci.values()),
            n_all_validated_platforms_both_supported=sum(
                all(any(int(z["both_supported"]) == 1 for z in entries if z["platform"] == p)
                    for p in set(z["platform"] for z in entries))
                for entries in loci.values()),
            all_manifest_runs_validated=int(len(complete_runs) == len(sruns))))
    sfields = ["sample", "species", "n_manifest_runs", "n_validated_runs",
               "platforms", "validated_platforms", "n_loci_with_results",
               "n_any_nuclear_supported", "n_any_chrM_supported", "n_any_both_supported",
               "n_all_validated_platforms_both_supported", "all_manifest_runs_validated"]
    write(outdir / "sample_species_summary.tsv", samples, sfields)
    print(f"Pilot: {len(samples)} samples / {len(set(x['species'] for x in samples))} species / {len(run_rows)} runs")
    print("Run statuses:", dict(Counter(x["status"] for x in run_rows)))
    print("Output:", outdir)


if __name__ == "__main__":
    main()
