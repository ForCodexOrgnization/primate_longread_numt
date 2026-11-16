#!/usr/bin/env python3
"""Independent post-alignment audit of chrM-homologous nuclear candidates.

No FASTQ download, no realignment. Inputs: indexed nuclear+chrM reference BAM,
reference FASTA, and besthit TSV. This measures sequence concordance *to the
nuclear reference*, not independent mtDNA origin or insertion discovery.
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path
import pysam


def numeric(value):
    try:
        return float(value)
    except (ValueError, TypeError):
        return None


def window_stats(pairs, sequence, reference, lo, hi):
    """Reference bases with aligned query bases, and exact base concordance."""
    aligned = matches = 0
    for qp, rp in pairs:
        if rp is None or rp < lo or rp >= hi or qp is None:
            continue
        aligned += 1
        if qp < len(sequence) and sequence[qp].upper() == reference[rp - lo].upper():
            matches += 1
    return aligned, matches


def maximal_deletion(read, lo, hi):
    """Largest D or N CIGAR event intersecting the nuclear candidate."""
    refpos = read.reference_start
    largest = 0
    for op, n in read.cigartuples or []:
        if op in (2, 3):
            overlap = min(refpos + n, hi) - max(refpos, lo)
            largest = max(largest, overlap)
        if op in (0, 2, 3, 7, 8):
            refpos += n
    return largest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bam", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--besthit", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--flank", type=int, default=100)
    p.add_argument("--min-flank-fraction", type=float, default=0.8)
    p.add_argument("--min-locus-fraction", type=float, default=0.9)
    p.add_argument("--min-identity", type=float, default=0.80)
    p.add_argument("--max-deletion", type=int, default=20)
    p.add_argument("--min-support", type=int, default=2)
    p.add_argument("--blast-identity", type=float, default=80)
    p.add_argument("--blast-length", type=int, default=100)
    args = p.parse_args()
    if args.flank < 1 or args.min_support < 1 or args.max_deletion < 0:
        p.error("flank/min-support must be positive and max-deletion nonnegative")
    if not all(0 < v <= 1 for v in (args.min_flank_fraction, args.min_locus_fraction, args.min_identity)):
        p.error("fraction and identity thresholds must lie in (0,1]")
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    with open(args.besthit, newline="") as fh:
        candidates = list(csv.DictReader(fh, delimiter="\t"))
    with pysam.FastaFile(args.reference) as fa, pysam.AlignmentFile(args.bam, "rb", reference_filename=args.reference) as bam:
        summary = []
        readrows = []
        for row in candidates:
            hit = str(row.get("chrM_hit", "")).strip()
            if not hit or hit in (".", "NA"):
                continue
            chrom = row["nuclear_chrom"]
            start, end = int(row["nuclear_start"]), int(row["nuclear_end"])
            if chrom not in fa.references or not (0 <= start < end <= fa.get_reference_length(chrom)):
                raise ValueError("Invalid candidate reference coordinates: " + str(row))
            identity = numeric(row.get("pident"))
            aln_length = numeric(row.get("align_length"))
            passing = (identity is not None and aln_length is not None and
                       identity >= args.blast_identity and aln_length >= args.blast_length)
            group = "BLAST_QC_PASS" if passing else "BLAST_QC_FAIL"
            left_lo = max(0, start - args.flank)
            right_hi = min(fa.get_reference_length(chrom), end + args.flank)
            left_ref = fa.fetch(chrom, left_lo, start)
            target_ref = fa.fetch(chrom, start, end)
            right_ref = fa.fetch(chrom, end, right_hi)
            counts = defaultdict(set)
            candidate_reads = 0
            for read in bam.fetch(chrom, left_lo, right_hi):
                if (read.is_unmapped or read.is_secondary or read.is_supplementary or
                    read.is_duplicate or read.is_qcfail or read.mapping_quality < args.min_mapq):
                    continue
                if read.query_sequence is None:
                    continue
                candidate_reads += 1
                pairs = read.get_aligned_pairs(matches_only=False)
                seq = read.query_sequence
                lt, lm = window_stats(pairs, seq, left_ref, left_lo, start)
                nt, nm = window_stats(pairs, seq, target_ref, start, end)
                rt, rm = window_stats(pairs, seq, right_ref, end, right_hi)
                left_ok = len(left_ref) == args.flank and lt >= args.min_flank_fraction * args.flank
                right_ok = len(right_ref) == args.flank and rt >= args.min_flank_fraction * args.flank
                fraction = nt / (end - start)
                seq_identity = nm / nt if nt else 0.0
                max_del = maximal_deletion(read, start, end)
                locus_ok = fraction >= args.min_locus_fraction and max_del <= args.max_deletion
                seq_ok = locus_ok and seq_identity >= args.min_identity
                both_ok = left_ok and right_ok and locus_ok
                strict_ok = both_ok and seq_ok
                name = read.query_name
                for flag, passed in [("left", left_ok), ("right", right_ok),
                                     ("locus", locus_ok), ("sequence", seq_ok),
                                     ("both", both_ok), ("strict", strict_ok)]:
                    if passed:
                        counts[flag].add(name)
                if left_ok or right_ok or nt > 0:
                    readrows.append({
                        "locus_id": row.get("locus_id", ""),
                        "group": group, "read": name,
                        "mapq": read.mapping_quality,
                        "left_aligned": lt, "locus_aligned": nt,
                        "right_aligned": rt, "locus_fraction": round(fraction, 5),
                        "locus_identity_to_nuclear_ref": round(seq_identity, 5),
                        "max_internal_deletion": max_del,
                        "left_ok": int(left_ok), "right_ok": int(right_ok),
                        "both_junction_windows_ok": int(both_ok),
                        "strict_sequence_support": int(strict_ok)
                    })
            summary.append({
                "locus_id": row.get("locus_id", ""),
                "nuclear_chrom": chrom, "nuclear_start0": start, "nuclear_end0": end,
                "chrM_hit": hit, "chrM_start": row.get("chrM_hit_start", row.get("chrM_start", "")),
                "chrM_end": row.get("chrM_hit_end", row.get("chrM_end", "")),
                "blast_pident": row.get("pident", ""),
                "blast_align_length": row.get("align_length", ""),
                "group": group, "n_primary_reads": candidate_reads,
                "n_left_window": len(counts["left"]),
                "n_right_window": len(counts["right"]),
                "n_continuous_locus": len(counts["locus"]),
                "n_sequence_concordant": len(counts["sequence"]),
                "n_both_junction_windows": len(counts["both"]),
                "n_strict_sequence_support": len(counts["strict"]),
                "strict_status": "SUPPORTED" if len(counts["strict"]) >= args.min_support else "NOT_SUPPORTED",
            })
    def write_tsv(path, rows, fields):
        with open(path, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
            wr.writeheader()
            wr.writerows(rows)
    locus_fields = ["locus_id", "nuclear_chrom", "nuclear_start0", "nuclear_end0",
                    "chrM_hit", "chrM_start", "chrM_end", "blast_pident",
                    "blast_align_length", "group", "n_primary_reads",
                    "n_left_window", "n_right_window", "n_continuous_locus",
                    "n_sequence_concordant", "n_both_junction_windows",
                    "n_strict_sequence_support", "strict_status"]
    read_fields = ["locus_id", "group", "read", "mapq", "left_aligned",
                   "locus_aligned", "right_aligned", "locus_fraction",
                   "locus_identity_to_nuclear_ref", "max_internal_deletion",
                   "left_ok", "right_ok", "both_junction_windows_ok",
                   "strict_sequence_support"]
    write_tsv(out / "chrM_candidate_locus_validation.tsv", summary, locus_fields)
    write_tsv(out / "chrM_candidate_read_evidence.tsv", readrows, read_fields)
    print("Candidates with recorded chrM hits:", len(summary))
    for group in ("BLAST_QC_PASS", "BLAST_QC_FAIL"):
        subset = [r for r in summary if r["group"] == group]
        supported = sum(r["strict_status"] == "SUPPORTED" for r in subset)
        print(f"{group}: {len(subset)} candidates, {supported} with >= {args.min_support} strict supporting reads")
    print("Outputs:", out)


if __name__ == "__main__":
    main()
