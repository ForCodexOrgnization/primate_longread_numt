#!/usr/bin/env python3
"""Plot strict nuclear NUMT support for ALL loci in completed pilot runs.

From repository root:
    python3 scripts/strict_nuclear_support_qc.py

Reads existing chrm_candidate_validation TSV files; does not touch BAMs or
change validation status. Both BLAST_QC_PASS and BLAST_QC_FAIL are included.
No-read loci are retained with zero support. Per-read MAPQ/identity panels are
conditional on read evidence: the upstream validator already excludes MAPQ<20,
secondary/supplementary, duplicates and QC-failed reads.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ORDER = ["Pongo_abelii", "Hoolock_leuconedys", "Colobus_guereza",
         "Rhinopithecus_roxellana", "Pygathrix_nigripes",
         "Ateles_geoffroyi", "Galago_moholi", "Loris_tardigradus"]
NUMERIC = ["mapq", "locus_identity_to_nuclear_ref", "locus_fraction",
           "left_aligned", "right_aligned", "max_internal_deletion",
           "strict_sequence_support"]
CRITERIA = [
    ("mapq", "MAPQ >=20*"),
    ("identity", "Identity >=80%"),
    ("locus_fraction", "Locus >=90%"),
    ("left", "Left flank >=80bp"),
    ("right", "Right flank >=80bp"),
    ("both_flanks", "Both flanks"),
    ("deletion", "Deletion <=20bp"),
    ("strict_2", ">=2 strict reads"),
]
METRICS = [
    ("best_mapq", "Best read MAPQ"),
    ("best_identity_pct", "Best read nuclear identity (%)"),
    ("best_locus_fraction_pct", "Best read locus fraction (%)"),
    ("best_left_bp", "Best read left aligned (bp)"),
    ("best_right_bp", "Best read right aligned (bp)"),
    ("n_strict_support_reads", "Strict supporting reads per locus"),
]


def write_plot(fig, out):
    fig.tight_layout()
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", default="results")
    p.add_argument("--summary", default="results/representative_pilot_summary")
    p.add_argument("--outdir", default=None)
    args = p.parse_args()
    root, summary = Path(args.root), Path(args.summary)
    out = Path(args.outdir) if args.outdir else summary / "strict_nuclear_support_qc"
    out.mkdir(parents=True, exist_ok=True)
    runs = pd.read_csv(summary / "run_summary.tsv", sep="\t", dtype=str).fillna("")
    runs = runs[runs["status"].eq("VALIDATED")]
    if runs.empty:
        raise SystemExit("No VALIDATED runs in run_summary.tsv")
    all_loci, all_reads = [], []
    for _, run in runs.iterrows():
        sample, accession = run["sample"], run["run"]
        species, platform = run["species"], run["platform"]
        d = root / sample / accession / "chrm_candidate_validation"
        lp, rp = d / "chrM_candidate_locus_validation.tsv", d / "chrM_candidate_read_evidence.tsv"
        if not (lp.is_file() and rp.is_file()):
            print("WARNING: missing evidence:", species, accession, flush=True)
            continue
        loci = pd.read_csv(lp, sep="\t", dtype={"locus_id": str}).fillna("")
        reads = pd.read_csv(rp, sep="\t", dtype={"locus_id": str}).fillna("")
        if loci["locus_id"].duplicated().any():
            raise ValueError(f"Duplicate locus IDs: {lp}")
        missing = set(NUMERIC) - set(reads.columns)
        if missing:
            raise ValueError(f"Missing read columns in {rp}: {sorted(missing)}")
        for field in NUMERIC:
            reads[field] = pd.to_numeric(reads[field], errors="coerce")
        reads["pass_mapq"] = reads["mapq"].ge(20)
        reads["pass_identity"] = reads["locus_identity_to_nuclear_ref"].ge(0.80)
        reads["pass_locus_fraction"] = reads["locus_fraction"].ge(0.90)
        reads["pass_left"] = reads["left_aligned"].ge(80)
        reads["pass_right"] = reads["right_aligned"].ge(80)
        reads["pass_both_flanks"] = reads["pass_left"] & reads["pass_right"]
        reads["pass_deletion"] = reads["max_internal_deletion"].le(20)
        reads["strict"] = reads["strict_sequence_support"].eq(1)
        for k, v in (("sample", sample), ("run", accession),
                     ("species", species), ("platform", platform)):
            reads[k] = v
        all_reads.append(reads)
        grouped = {str(k): g for k, g in reads.groupby("locus_id", sort=False)}
        for _, loc in loci.iterrows():
            lid = str(loc["locus_id"])
            rr = grouped.get(lid)
            best = None
            if rr is not None and not rr.empty:
                best = rr.sort_values(["locus_fraction", "locus_identity_to_nuclear_ref", "mapq"],
                                      ascending=False, na_position="last").iloc[0]
            def bestval(field, scale=1):
                return float(best[field]) * scale if best is not None and pd.notna(best[field]) else np.nan
            def anypass(field):
                return bool(rr["pass_" + field].any()) if rr is not None and not rr.empty else False
            n = int(pd.to_numeric(loc["n_strict_sequence_support"], errors="raise"))
            item = dict(sample=sample, run=accession, species=species, platform=platform,
                        locus_id=lid, group=loc["group"], n_primary_reads=int(loc["n_primary_reads"]),
                        n_reads_recorded=(rr["read"].nunique() if rr is not None else 0),
                        n_strict_support_reads=n, strict_supported_locus=n >= 2,
                        best_mapq=bestval("mapq"),
                        best_identity_pct=bestval("locus_identity_to_nuclear_ref", 100),
                        best_locus_fraction_pct=bestval("locus_fraction", 100),
                        best_left_bp=bestval("left_aligned"),
                        best_right_bp=bestval("right_aligned"),
                        best_max_internal_deletion=bestval("max_internal_deletion"))
            for field in ("mapq", "identity", "locus_fraction", "left",
                          "right", "both_flanks", "deletion"):
                item["any_" + field] = anypass(field)
            all_loci.append(item)
        print(f"{species}: {len(loci)} loci, {len(reads)} read-evidence rows", flush=True)
    if not all_loci:
        raise SystemExit("No per-locus results collected")
    L = pd.DataFrame(all_loci)
    R = pd.concat(all_reads, ignore_index=True)
    L.to_csv(out / "all_locus_summary.tsv", sep="\t", index=False)
    R.to_csv(out / "all_read_evidence_combined.tsv", sep="\t", index=False)

    summaries = []
    for (species, platform), g in L.groupby(["species", "platform"], sort=False):
        item = dict(species=species, platform=platform, n_loci=len(g),
                    n_blast_pass=int(g["group"].eq("BLAST_QC_PASS").sum()),
                    n_strict_supported=int(g["strict_supported_locus"].sum()),
                    median_strict_reads=float(g["n_strict_support_reads"].median()))
        for field, _ in CRITERIA:
            passed = g["strict_supported_locus"] if field == "strict_2" else g["any_" + field]
            item["n_" + field] = int(passed.sum())
            item["pct_" + field] = round(100 * passed.mean(), 2)
        summaries.append(item)
    C = pd.DataFrame(summaries)
    C["_order"] = C["species"].map({s: i for i, s in enumerate(ORDER)}).fillna(999)
    C = C.sort_values("_order").drop(columns="_order")
    C.to_csv(out / "criterion_pass_summary_by_species.tsv", sep="\t", index=False)

    order = C["species"].tolist()
    for field, ylabel in METRICS:
        values = [L.loc[L["species"].eq(s), field].dropna().to_numpy(float) for s in order]
        fig, ax = plt.subplots(figsize=(12, 5))
        ax.boxplot([v if len(v) else np.array([np.nan]) for v in values], showfliers=False)
        ax.set_xticks(range(1, len(order) + 1), [s.replace("_", " ") for s in order], rotation=35, ha="right")
        ax.set_ylabel(ylabel)
        ax.set_title("Per-locus distribution (all candidate NUMTs)")
        ax.grid(axis="y", alpha=.2)
        write_plot(fig, out / (field + ".boxplot.png"))

    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(C))
    width = .095
    for i, (field, label) in enumerate(CRITERIA):
        ax.bar(x + (i - (len(CRITERIA)-1)/2)*width, C["pct_" + field], width, label=label)
    ax.set_xticks(x, [s.replace("_", " ") for s in order], rotation=35, ha="right")
    ax.set_ylim(0, 105)
    ax.set_ylabel("% of candidate NUMT loci")
    ax.set_title("Independent criterion pass rates (not a sequential funnel)")
    ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(.5,-.22), frameon=False)
    write_plot(fig, out / "criterion_pass_rates_by_species.png")

    labels = ["0", "1", "2", "3-5", ">5"]
    L["support_bin"] = pd.cut(L["n_strict_support_reads"],
                              bins=[-.1, .5, 1.5, 2.5, 5.5, np.inf], labels=labels)
    T = pd.crosstab(L["species"], L["support_bin"]).reindex(index=order, columns=labels, fill_value=0)
    T.div(T.sum(axis=1), axis=0).mul(100).to_csv(
        out / "strict_support_read_bins_by_species.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(12, 5))
    bottom = np.zeros(len(order))
    for field in labels:
        vals = T[field].to_numpy() / T.sum(axis=1).to_numpy() * 100
        ax.bar(np.arange(len(order)), vals, bottom=bottom, label=field)
        bottom += vals
    ax.set_xticks(np.arange(len(order)), [s.replace("_", " ") for s in order], rotation=35, ha="right")
    ax.set_ylim(0, 100)
    ax.set_ylabel("% of candidate NUMT loci")
    ax.set_title("Strict nuclear-supporting reads per NUMT")
    ax.legend(title="# reads", ncol=5, frameon=False)
    write_plot(fig, out / "strict_support_read_count_distribution.png")
    print("\n", C[["species", "n_loci", "n_blast_pass", "n_strict_supported", "median_strict_reads"]].to_string(index=False))
    print("\nOutputs:", out)
    print("NOTE: MAPQ rates are conditional on upstream MAPQ>=20 filtering; this script cannot test lower MAPQ thresholds.")


if __name__ == "__main__":
    main()
