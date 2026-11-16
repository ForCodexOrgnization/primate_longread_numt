#!/usr/bin/env python3
"""BLAST nuclear-mapped ONT read segments against species-matched chrM.

Requires existing validate_chrm_candidates.py output (read evidence + locus
validation), so only reads passing the same strict nuclear-locus checks are
tested. This is a targeted homology check, not de novo NUMT assembly.
"""
import argparse
import csv
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import pysam


def read_tsv(path):
    with open(path, newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def integer(value):
    try:
        return int(float(value))
    except (ValueError, TypeError):
        return None


def extract_aligned_segment(read, start, end):
    """Read bases aligned in nuclear [start,end), including internal insertions.

    pysam read.query_sequence is in alignment orientation for reverse alignments.
    Soft clips are excluded; reference deletions contribute no query bases.
    """
    reference_pos = read.reference_start
    query_pos = 0
    pieces = []
    for op, size in read.cigartuples or []:
        if op in (0, 7, 8):  # M, =, X
            lo = max(start, reference_pos)
            hi = min(end, reference_pos + size)
            if lo < hi:
                pieces.append(read.query_sequence[query_pos + lo - reference_pos:
                                                  query_pos + hi - reference_pos])
            reference_pos += size
            query_pos += size
        elif op == 1:  # inserted bases in query
            if start < reference_pos < end:
                pieces.append(read.query_sequence[query_pos:query_pos + size])
            query_pos += size
        elif op in (2, 3):  # D, N
            reference_pos += size
        elif op == 4:  # S (query only)
            query_pos += size
        elif op in (5, 6):  # H, P
            pass
        else:
            raise ValueError(f"Unexpected CIGAR opcode {op}")
    return "".join(pieces).upper()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", required=True)
    p.add_argument("--reference", required=True, help="FASTA containing chrM")
    p.add_argument("--nuclear-validation", required=True, help="Existing chrM_candidate_validation directory")
    p.add_argument("--outdir", required=True)
    p.add_argument("--chrm-contig", default="chrM")
    p.add_argument("--max-reads-per-locus", type=int, default=40)
    p.add_argument("--min-blast-identity", type=float, default=75.0)
    p.add_argument("--min-blast-length", type=int, default=75)
    p.add_argument("--expected-margin", type=int, default=250)
    p.add_argument("--min-support", type=int, default=2)
    p.add_argument("--threads", type=int, default=2)
    args = p.parse_args()
    if min(args.max_reads_per_locus, args.min_blast_length, args.min_support, args.threads) <= 0:
        p.error("Read limit, alignment length, support and threads must be positive")
    if not 0 < args.min_blast_identity <= 100 or args.expected_margin < 0:
        p.error("Invalid BLAST identity or margin")
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    validated = Path(args.nuclear_validation)
    loci = read_tsv(validated / "chrM_candidate_locus_validation.tsv")
    nuclear_reads = read_tsv(validated / "chrM_candidate_read_evidence.tsv")
    candidates = {}
    for row in loci:
        key = row["locus_id"]
        if key in candidates:
            raise ValueError("Duplicate locus_id: " + key)
        candidates[key] = row
    strict = defaultdict(dict)
    for row in nuclear_reads:
        if row["locus_id"] not in candidates or row["strict_sequence_support"] != "1":
            continue
        name = row["read"]
        key = row["locus_id"]
        # Read names are only unique within one locus here.
        if name not in strict[key] or int(row["mapq"]) > int(strict[key][name]["mapq"]):
            strict[key][name] = row

    selected = {}
    for key, reads in strict.items():
        ranked = sorted(reads.values(), key=lambda r: (
            -int(r["mapq"]), -float(r["locus_fraction"]), r["read"]))
        selected[key] = {r["read"] for r in ranked[:args.max_reads_per_locus]}

    query_map = {}
    extraction_counts = Counter()
    queries = out / "nuclear_mapped_read_segments.fa"
    with pysam.FastaFile(args.reference) as ref:
        if args.chrm_contig not in ref.references:
            raise ValueError(f"Reference lacks chrM contig: {args.chrm_contig}")
        chrm_fa = out / "chrM_reference.fa"
        chrm_fa.write_text(f">{args.chrm_contig}\n{ref.fetch(args.chrm_contig)}\n")
    with pysam.AlignmentFile(args.bam, "rb", reference_filename=args.reference) as bam, open(queries, "w") as fasta:
        for key, row in candidates.items():
            names = selected.get(key, set())
            if not names:
                continue
            chrom = row["nuclear_chrom"]
            start = int(row["nuclear_start0"])
            end = int(row["nuclear_end0"])
            found = set()
            for read in bam.fetch(chrom, start, end):
                if (read.query_name not in names or read.query_name in found
                    or read.is_secondary or read.is_supplementary
                    or read.is_unmapped or read.is_duplicate or read.is_qcfail
                    or read.query_sequence is None):
                    continue
                seq = extract_aligned_segment(read, start, end)
                if len(seq) < args.min_blast_length:
                    continue
                found.add(read.query_name)
                qid = f"q{len(query_map):07d}"
                query_map[qid] = (key, read.query_name, len(seq))
                fasta.write(f">{qid}\n{seq}\n")
                extraction_counts[key] += 1

    fields = ["query_id", "locus_id", "read", "group", "segment_length",
              "best_chrM_pident", "best_chrM_align_length", "best_chrM_start",
              "best_chrM_end", "any_homology", "expected_chrM_homology",
              "expected_pident", "expected_align_length", "expected_start",
              "expected_end"]
    raw_path = out / "read_segments_vs_chrM.blast.tsv"
    hit_fields = ["qseqid", "sseqid", "pident", "length", "qstart", "qend",
                  "sstart", "send", "bitscore", "evalue"]
    hits = defaultdict(list)
    if query_map:
        with open(raw_path, "w") as handle:
            subprocess.run([
                "blastn", "-task", "blastn", "-dust", "no", "-word_size", "7",
                "-query", str(queries), "-subject", str(chrm_fa),
                "-num_threads", str(args.threads), "-evalue", "1e-2",
                "-outfmt", "6 " + " ".join(hit_fields),
            ], stdout=handle, check=True)
        for h in read_tsv_blast(raw_path, hit_fields):
            if h["qseqid"] in query_map:
                hits[h["qseqid"]].append(h)
    else:
        raw_path.write_text("")

    read_out = []
    by_locus = defaultdict(lambda: {"any": set(), "expected": set()})
    for qid, (key, name, segment_len) in query_map.items():
        locus = candidates[key]
        all_hits = hits[qid]
        qualifying = [h for h in all_hits if float(h["pident"]) >= args.min_blast_identity
                      and int(h["length"]) >= args.min_blast_length]
        qualifying.sort(key=lambda h: float(h["bitscore"]), reverse=True)
        exp_start = integer(locus.get("chrM_start"))
        exp_end = integer(locus.get("chrM_end"))
        expected = []
        if exp_start is not None and exp_end is not None:
            left, right = sorted((exp_start, exp_end))
            left = max(1, left - args.expected_margin)
            right += args.expected_margin
            for h in qualifying:
                hs, he = sorted((int(h["sstart"]), int(h["send"])))
                if min(right, he) >= max(left, hs):
                    expected.append(h)
        best = qualifying[0] if qualifying else None
        best_expected = expected[0] if expected else None
        if best:
            by_locus[key]["any"].add(name)
        if best_expected:
            by_locus[key]["expected"].add(name)
        read_out.append({
            "query_id": qid, "locus_id": key, "read": name,
            "group": locus["group"], "segment_length": segment_len,
            "best_chrM_pident": best["pident"] if best else "",
            "best_chrM_align_length": best["length"] if best else "",
            "best_chrM_start": min(int(best["sstart"]), int(best["send"])) if best else "",
            "best_chrM_end": max(int(best["sstart"]), int(best["send"])) if best else "",
            "any_homology": int(bool(best)),
            "expected_chrM_homology": int(bool(best_expected)),
            "expected_pident": best_expected["pident"] if best_expected else "",
            "expected_align_length": best_expected["length"] if best_expected else "",
            "expected_start": min(int(best_expected["sstart"]), int(best_expected["send"])) if best_expected else "",
            "expected_end": max(int(best_expected["sstart"]), int(best_expected["send"])) if best_expected else "",
        })
    write_tsv(out / "per_read_chrM_homology.tsv", read_out, fields)
    locus_fields = ["locus_id", "group", "n_nuclear_strict_reads",
                    "n_tested_read_segments", "n_any_chrM_homology",
                    "n_expected_chrM_homology", "chrM_status",
                    "chrM_hit", "expected_chrM_start", "expected_chrM_end",
                    "blast_pident", "blast_align_length"]
    locus_out = []
    for key, locus in candidates.items():
        n_any = len(by_locus[key]["any"])
        n_expected = len(by_locus[key]["expected"])
        locus_out.append({
            "locus_id": key, "group": locus["group"],
            "n_nuclear_strict_reads": locus["n_strict_sequence_support"],
            "n_tested_read_segments": extraction_counts[key],
            "n_any_chrM_homology": n_any,
            "n_expected_chrM_homology": n_expected,
            "chrM_status": ("EXPECTED_REGION_SUPPORTED" if n_expected >= args.min_support
                            else "OTHER_REGION_ONLY" if n_any >= args.min_support
                            else "NOT_SUPPORTED"),
            "chrM_hit": locus["chrM_hit"],
            "expected_chrM_start": locus["chrM_start"],
            "expected_chrM_end": locus["chrM_end"],
            "blast_pident": locus["blast_pident"],
            "blast_align_length": locus["blast_align_length"],
        })
    write_tsv(out / "per_locus_chrM_homology.tsv", locus_out, locus_fields)
    print(f"[chrM homology] tested {len(query_map)} unique locus/read sequences across {len(candidates)} loci")
    for group in ("BLAST_QC_PASS", "BLAST_QC_FAIL"):
        rows = [r for r in locus_out if r["group"] == group]
        counts = Counter(r["chrM_status"] for r in rows)
        print(f"{group}: N={len(rows)}; " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    print("Outputs:", out)


def read_tsv_blast(path, fields):
    with open(path) as stream:
        for line in stream:
            if not line.strip():
                continue
            values = line.rstrip("\n").split("\t")
            if len(values) != len(fields):
                raise ValueError("Unexpected BLAST output row: " + line)
            yield dict(zip(fields, values))


if __name__ == "__main__":
    main()
