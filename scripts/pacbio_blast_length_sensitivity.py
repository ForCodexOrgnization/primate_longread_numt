#!/usr/bin/env python3
"""Exploratory PacBio NUMT chrM BLAST-length sensitivity analysis.

Run at repository root:
  python3 scripts/pacbio_blast_length_sensitivity.py

Requires: pandas, pysam, BLAST+ (blastn, makeblastdb), and existing completed
PacBio BAM/nuclear-validation/chrM-homology outputs. Original outputs are never
modified. This repeats BLAST with a 30-bp extraction floor, because the
original pipeline never tested segments shorter than 75 bp.
"""
import argparse
import math
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
import pysam

THRESHOLDS = (30, 50, 75, 100)
BINS = ((0, 100, "<100"), (100, 200, "100-199"),
        (200, 500, "200-499"), (500, 1000, "500-999"),
        (1000, float("inf"), ">=1000"))


def extract_segment(read, start, end):
    """Extract query bases aligned to nuclear interval, including insertions."""
    refpos = read.reference_start
    querypos = 0
    pieces = []
    for op, n in read.cigartuples or []:
        if op in (0, 7, 8):
            lo, hi = max(start, refpos), min(end, refpos + n)
            if lo < hi:
                pieces.append(read.query_sequence[querypos + lo - refpos:
                                                   querypos + hi - refpos])
            refpos += n
            querypos += n
        elif op == 1:
            if start < refpos < end:
                pieces.append(read.query_sequence[querypos:querypos + n])
            querypos += n
        elif op in (2, 3):
            refpos += n
        elif op == 4:
            querypos += n
        elif op not in (5, 6):
            raise ValueError(f"Unexpected CIGAR op {op}")
    return "".join(pieces).upper()


def numeric(value, default=None):
    try:
        x = float(value)
        return x if math.isfinite(x) else default
    except (TypeError, ValueError):
        return default


def overlap(a, b, c, d):
    return max(0, min(b, d) - max(a, c))


def length_bin(length):
    return next((label for lo, hi, label in BINS if lo <= length < hi), "unknown")


def one_run(run, root, args, temp):
    sample, accession, species = map(str, (run["sample"], run["run"], run["species"]))
    base = root / sample / accession
    nuclear_dir = base / "chrm_candidate_validation"
    homology_dir = base / "longread_chrM_homology"
    bamfile = base / "aligned.bam"
    loci_file = nuclear_dir / "chrM_candidate_locus_validation.tsv"
    reads_file = nuclear_dir / "chrM_candidate_read_evidence.tsv"
    chrm_file = homology_dir / "chrM_reference.fa"
    missing = [str(p) for p in (bamfile, loci_file, reads_file, chrm_file)
               if not p.is_file()]
    if missing:
        print(f"SKIP {species} {accession}: missing {missing}", flush=True)
        return [], [], []
    loci = pd.read_csv(loci_file, sep="\t", dtype={"locus_id": str}).fillna("")
    reads = pd.read_csv(reads_file, sep="\t", dtype={"locus_id": str}).fillna("")
    if args.group != "ALL":
        loci = loci[loci["group"] == args.group].copy()
    if loci.empty:
        return [], [], []
    strict = reads[pd.to_numeric(reads["strict_sequence_support"], errors="coerce") == 1].copy()
    strict["mapq_sort"] = pd.to_numeric(strict["mapq"], errors="coerce").fillna(-1)
    strict["frac_sort"] = pd.to_numeric(strict["locus_fraction"], errors="coerce").fillna(0)
    strict = strict.sort_values(["mapq_sort", "frac_sort"], ascending=False)
    strict = strict.drop_duplicates(["locus_id", "read"])
    by_locus = {str(key): set(group.head(args.max_reads_per_locus)["read"].astype(str))
                for key, group in strict.groupby("locus_id")}
    expected = {}
    for _, r in loci.iterrows():
        a, b = numeric(r.get("chrM_start")), numeric(r.get("chrM_end"))
        expected[str(r["locus_id"])] = ((max(0, int(min(a, b)) - 1), int(max(a, b)))
                                          if a is not None and b is not None else None)

    qmap = {}
    queryfile = temp / f"{accession}.queries.fa"
    with pysam.AlignmentFile(str(bamfile), "rb") as bam, queryfile.open("w") as handle:
        for _, loc in loci.iterrows():
            lid = str(loc["locus_id"])
            names = by_locus.get(lid, set())
            if not names:
                continue
            chrom = str(loc["nuclear_chrom"])
            start, end = int(loc["nuclear_start0"]), int(loc["nuclear_end0"])
            found = set()
            for read in bam.fetch(chrom, start, end):
                if (read.query_name not in names or read.query_name in found
                    or read.is_unmapped or read.is_secondary or read.is_supplementary
                    or read.is_duplicate or read.is_qcfail or read.query_sequence is None):
                    continue
                segment = extract_segment(read, start, end)
                found.add(read.query_name)
                if len(segment) < min(THRESHOLDS):
                    continue
                qid = f"q{len(qmap):09d}"
                qmap[qid] = (lid, read.query_name, len(segment))
                handle.write(f">{qid}\n{segment}\n")

    hits = defaultdict(list)
    if qmap:
        database = temp / f"{accession}.chrMdb"
        subprocess.run(["makeblastdb", "-in", str(chrm_file), "-dbtype", "nucl",
                        "-out", str(database)], check=True, stdout=subprocess.DEVNULL)
        output = temp / f"{accession}.blast.tsv"
        with output.open("w") as fh:
            subprocess.run(["blastn", "-task", "blastn", "-dust", "no",
                            "-word_size", str(args.word_size), "-query", str(queryfile),
                            "-db", str(database), "-num_threads", str(args.threads),
                            "-evalue", str(args.evalue),
                            "-outfmt", "6 qseqid pident length sstart send bitscore"],
                           stdout=fh, check=True)
        with output.open() as fh:
            for line in fh:
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 6 or fields[0] not in qmap:
                    continue
                qid, pct, alnlen, start, end, score = fields
                hits[qid].append((float(pct), int(alnlen),
                                  min(int(start), int(end))-1,
                                  max(int(start), int(end)), float(score)))

    q_by_locus = defaultdict(list)
    for qid, (lid, name, seglen) in qmap.items():
        q_by_locus[lid].append((qid, name, seglen))
    hit_rows, locus_rows, threshold_rows = [], [], []
    for _, loc in loci.iterrows():
        lid = str(loc["locus_id"])
        start, end = int(loc["nuclear_start0"]), int(loc["nuclear_end0"])
        iv = expected[lid]
        rec = dict(sample=sample, run=accession, species=species,
                   locus_id=lid, group=str(loc["group"]),
                   numt_length=end-start, length_bin=length_bin(end-start),
                   n_nuclear_strict=int(numeric(loc["n_strict_sequence_support"], 0)),
                   n_tested_segments_30bp=len(q_by_locus[lid]),
                   original_chrM_start1=iv[0]+1 if iv else "",
                   original_chrM_end1=iv[1] if iv else "")
        locus_rows.append(rec)
        for qid, name, seglen in q_by_locus[lid]:
            for pct, alen, hs, he, score in hits[qid]:
                orig = overlap(hs, he, *iv) if iv else 0
                expanded = (max(0, iv[0]-args.margin), iv[1]+args.margin) if iv else None
                margin = overlap(hs, he, *expanded) if expanded else 0
                hit_rows.append({**rec, "read":name, "segment_length":seglen,
                                 "blast_identity_pct":pct, "blast_align_length":alen,
                                 "blast_chrM_start1":hs+1, "blast_chrM_end1":he,
                                 "blast_bitscore":score, "overlap_original_bp":orig,
                                 "overlap_expanded_bp":margin,
                                 "fraction_hit_overlapping_original":round(orig/(he-hs), 5),
                                 "fraction_original_covered":round(orig/(iv[1]-iv[0]), 5) if iv else "",
                                 "overlap_class":"ORIGINAL" if orig else
                                                 "MARGIN_ONLY" if margin else "OTHER_REGION"})
        for minimum in THRESHOLDS:
            groups = {"ORIGINAL":set(), "MARGIN_ONLY":set(), "OTHER_REGION":set()}
            qualifying_hsps = 0
            for qid, name, _ in q_by_locus[lid]:
                valid = [h for h in hits[qid] if h[0] >= args.identity and h[1] >= minimum]
                qualifying_hsps += len(valid)
                if not valid:
                    continue
                if iv and any(overlap(h[2],h[3],*iv)>0 for h in valid):
                    groups["ORIGINAL"].add(name)
                elif iv and any(overlap(h[2],h[3],max(0,iv[0]-args.margin),iv[1]+args.margin)>0 for h in valid):
                    groups["MARGIN_ONLY"].add(name)
                else:
                    groups["OTHER_REGION"].add(name)
            n_orig = len(groups["ORIGINAL"])
            n_margin = len(groups["MARGIN_ONLY"])
            n_expanded = n_orig+n_margin
            threshold_rows.append({**rec, "threshold_bp":minimum,
                                   "identity_threshold_pct":args.identity,
                                   "n_expected_original_reads":n_orig,
                                   "n_margin_only_reads":n_margin,
                                   "n_expected_expanded_reads":n_expanded,
                                   "n_other_region_only_reads":len(groups["OTHER_REGION"]),
                                   "n_qualifying_hsps":qualifying_hsps,
                                   "pass_nuclear_2":int(rec["n_nuclear_strict"]>=2),
                                   "pass_original_chrM_2":int(n_orig>=2),
                                   "pass_expanded_chrM_2":int(n_expanded>=2),
                                   "pass_both_original":int(rec["n_nuclear_strict"]>=2 and n_orig>=2),
                                   "pass_both_expanded":int(rec["n_nuclear_strict"]>=2 and n_expanded>=2)})
    print(f"{species}: {len(loci)} loci; {len(qmap)} read segments; {len(hit_rows)} BLAST HSPs",
          flush=True)
    return hit_rows, locus_rows, threshold_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default="results")
    parser.add_argument("--summary", default="results/representative_pilot_summary")
    parser.add_argument("--identity", type=float, default=75)
    parser.add_argument("--margin", type=int, default=250)
    parser.add_argument("--evalue", type=float, default=1e-2)
    parser.add_argument("--word-size", type=int, default=7)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--max-reads-per-locus", type=int, default=40)
    parser.add_argument("--group", choices=["BLAST_QC_PASS","BLAST_QC_FAIL","ALL"],
                        default="BLAST_QC_PASS")
    args = parser.parse_args()
    root, summary = Path(args.root), Path(args.summary)
    runs = pd.read_csv(summary/"run_summary.tsv", sep="\t").fillna("")
    runs = runs[(runs["status"]=="VALIDATED") &
                (runs["platform"].astype(str).str.upper()=="PACBIO")]
    if runs.empty:
        raise SystemExit("No validated PacBio runs in run_summary.tsv")
    all_hits, all_loci, all_thresholds = [], [], []
    with tempfile.TemporaryDirectory(prefix="pacbio_numt_blast_") as temp:
        for _, r in runs.iterrows():
            h, l, t = one_run(r, root, args, Path(temp))
            all_hits += h
            all_loci += l
            all_thresholds += t
    if not all_loci:
        raise SystemExit("No PacBio inputs found")
    output = summary/"pacbio_blast_length_sensitivity"
    output.mkdir(parents=True, exist_ok=True)
    hdf, ldf, tdf = (pd.DataFrame(all_hits), pd.DataFrame(all_loci),
                      pd.DataFrame(all_thresholds))
    hdf.to_csv(output/"all_blast_hsps.tsv", sep="\t", index=False)
    ldf.to_csv(output/"loci.tsv", sep="\t", index=False)
    tdf.to_csv(output/"per_locus_thresholds.tsv", sep="\t", index=False)
    aggregations = dict(n_loci=("locus_id","size"),
                        n_nuclear_ge2=("pass_nuclear_2","sum"),
                        n_original_ge2=("pass_original_chrM_2","sum"),
                        n_expanded_ge2=("pass_expanded_chrM_2","sum"),
                        n_both_original=("pass_both_original","sum"),
                        n_both_expanded=("pass_both_expanded","sum"))
    by_length = tdf.groupby(["species","length_bin","threshold_bp"],as_index=False).agg(**aggregations)
    by_length["validation_rate_expanded_pct"] = (100*by_length.n_both_expanded/by_length.n_loci).round(2)
    by_length.to_csv(output/"by_species_length_bin.tsv", sep="\t", index=False)
    by_species = tdf.groupby(["species","threshold_bp"],as_index=False).agg(**aggregations)
    by_species = by_species.rename(columns={"n_original_ge2":"n_chrM_original_ge2",
                                             "n_expanded_ge2":"n_chrM_expanded_ge2"})
    by_species["validation_rate_expanded_pct"] = (100*by_species.n_both_expanded/by_species.n_loci).round(2)
    by_species.to_csv(output/"by_species_threshold.tsv", sep="\t", index=False)
    print("\n=== PacBio chrM BLAST length sensitivity ===")
    print(by_species.to_string(index=False))
    if not hdf.empty:
        print("\n=== BLAST HSP lengths by species and region match ===")
        print(hdf.groupby(["species","overlap_class"])["blast_align_length"].agg(["count","median","min","max"]).to_string())
    print("\nOutputs:", output)
    print("Exploratory: strict nuclear read selection is unchanged; 30-bp alignments may be nonspecific.")


if __name__ == "__main__":
    main()
