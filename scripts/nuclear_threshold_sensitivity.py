#!/usr/bin/env python3
"""Read-only nuclear NUMT threshold sensitivity analysis across completed runs.

Run from repo root:
  python3 scripts/nuclear_threshold_sensitivity.py

Uses existing per-read evidence; no BAM realignment or original QC changes.
Includes ALL BLAST_QC_PASS and BLAST_QC_FAIL NUMT loci, including zero-read loci.
MAPQ<20 reads and reads absent from the evidence table CANNOT be recovered.
This is a sensitivity analysis, not a new validated callset.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ORDER = ["Pongo_abelii","Hoolock_leuconedys","Colobus_guereza",
         "Rhinopithecus_roxellana","Pygathrix_nigripes",
         "Ateles_geoffroyi","Galago_moholi","Loris_tardigradus"]
REQUIRED = ["locus_id","read","mapq","left_aligned","right_aligned",
            "locus_fraction","locus_identity_to_nuclear_ref","max_internal_deletion"]
NUCLEAR_COLS = ["nuclear_start0","nuclear_end0"]


def scenarios():
    original = dict(flank=80, fraction=.90, reads=2, deletion=20)
    configs = [("Original", original),
               ("Moderate", dict(flank=50,fraction=.80,reads=2,deletion=20)),
               ("Exploratory", dict(flank=50,fraction=.80,reads=1,deletion=20))]
    for n in (50,30):
        configs.append((f"Flank_{n}",{**original,"flank":n}))
    for value in (.80,.70):
        configs.append((f"Locus_{round(value*100)}",{**original,"fraction":value}))
    configs.append(("Reads_1",{**original,"reads":1}))
    configs.append(("Reads_3",{**original,"reads":3}))
    configs.append(("Deletion_50",{**original,"deletion":50}))
    return configs


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root",default="results")
    p.add_argument("--summary",default="results/representative_pilot_summary")
    p.add_argument("--outdir",default=None)
    args=p.parse_args()
    root=Path(args.root)
    summary=Path(args.summary)
    out=Path(args.outdir) if args.outdir else summary/"nuclear_threshold_sensitivity"
    out.mkdir(parents=True,exist_ok=True)
    runs=pd.read_csv(summary/"run_summary.tsv",sep="\t",dtype=str).fillna("")
    runs=runs[runs.status.eq("VALIDATED")].copy()
    if runs.empty:
        raise SystemExit("No VALIDATED runs found")

    all_loci=[]
    all_details=[]
    configs=scenarios()
    for _,run in runs.iterrows():
        sample, accession, species, platform=[str(run[k]) for k in ("sample","run","species","platform")]
        d=root/sample/accession/"chrm_candidate_validation"
        lp=d/"chrM_candidate_locus_validation.tsv"
        rp=d/"chrM_candidate_read_evidence.tsv"
        if not (lp.exists() and rp.exists()):
            print(f"WARNING missing evidence for {species}: {d}",flush=True)
            continue
        loc=pd.read_csv(lp,sep="\t",dtype={"locus_id":str}).fillna("")
        read=pd.read_csv(rp,sep="\t",dtype={"locus_id":str}).fillna("")
        missing=set(REQUIRED)-set(read.columns)
        if missing:
            raise ValueError(f"Missing fields in {rp}: {sorted(missing)}")
        if loc.locus_id.duplicated().any():
            raise ValueError(f"Duplicate loci in {lp}")
        for col in REQUIRED[2:]:
            read[col]=pd.to_numeric(read[col],errors="coerce")
        # Preserve original validator requirement: full 100 bp windows on both sides.
        # The saved per-read evidence does not explicitly record flank-window lengths.
        # Loci adjacent to contig boundaries cannot be assessed reliably from TSV only.
        loc["window_assessable"] = True
        if set(NUCLEAR_COLS).issubset(loc.columns):
            pos=pd.to_numeric(loc["nuclear_start0"],errors="coerce")
            loc["window_assessable"] = pos.ge(100)
        # Count UNIQUE read names per locus (validator uses sets).
        groups={str(k):g for k,g in read.groupby("locus_id",sort=False)}
        for _,lr in loc.iterrows():
            lid=str(lr["locus_id"])
            g=groups.get(lid)
            if g is None:
                g=read.iloc[:0]
            baseline=int(float(lr["n_strict_sequence_support"]))
            for name,cfg in configs:
                mask=(g["mapq"].ge(20)
                      & g["locus_identity_to_nuclear_ref"].ge(.80)
                      & g["locus_fraction"].ge(cfg["fraction"])
                      & g["left_aligned"].ge(cfg["flank"])
                      & g["right_aligned"].ge(cfg["flank"])
                      & g["max_internal_deletion"].le(cfg["deletion"]))
                if not bool(lr["window_assessable"]):
                    mask=mask & False
                n=int(g.loc[mask,"read"].nunique())
                detail=dict(species=species,platform=platform,sample=sample,run=accession,
                            locus_id=lid,group=lr["group"],scenario=name,
                            flank_bp=cfg["flank"],locus_fraction=cfg["fraction"],
                            max_deletion_bp=cfg["deletion"],min_support_reads=cfg["reads"],
                            n_eligible_reads=n,pass_criterion=(n>=cfg["reads"]),
                            original_strict_reads=baseline,original_pass=(baseline>=2),
                            n_primary_reads=int(float(lr["n_primary_reads"])),
                            window_assessable=bool(lr["window_assessable"]))
                all_details.append(detail)
                if name=="Original" and n!=baseline:
                    print(f"WARNING original mismatch {species}/{lid}: recomputed={n} recorded={baseline}",flush=True)
            all_loci.append(dict(species=species,platform=platform,sample=sample,run=accession,
                                 locus_id=lid,group=lr["group"],n_primary_reads=int(float(lr["n_primary_reads"])),
                                 n_original_strict=baseline,original_pass=baseline>=2))
        print(f"{species} {accession}: {len(loc)} candidates",flush=True)
    if not all_details:
        raise SystemExit("No results generated")
    details=pd.DataFrame(all_details)
    details.to_csv(out/"per_locus_scenarios.tsv",sep="\t",index=False)
    S=(details.groupby(["species","platform","group","scenario"],as_index=False)
       .agg(n_loci=("locus_id","size"),n_supported=("pass_criterion","sum"),
            n_original_supported=("original_pass","sum")))
    S["rate_pct"]=(S.n_supported/S.n_loci*100).round(2)
    S["new_supported_vs_original"]=S.n_supported-S.n_original_supported
    S.to_csv(out/"by_species_group_scenario.tsv",sep="\t",index=False)
    A=(details.groupby(["species","platform","scenario"],as_index=False)
       .agg(n_loci=("locus_id","size"),n_supported=("pass_criterion","sum"),
            n_original_supported=("original_pass","sum")))
    A["rate_pct"]=(A.n_supported/A.n_loci*100).round(2)
    A["new_supported_vs_original"]=A.n_supported-A.n_original_supported
    A.to_csv(out/"by_species_all_candidates.tsv",sep="\t",index=False)

    # Focused Colobus report: include both BLAST QC groups, and flag changes per scenario.
    target=details[details.species.eq("Colobus_guereza")].copy()
    failed=target[target.scenario.eq("Original") & (~target.original_pass)]
    target=target[target.locus_id.isin(failed.locus_id)]
    target.to_csv(out/"colobus_original_failures_all_scenarios.tsv",sep="\t",index=False)
    print("\n=== ALL NUMT CANDIDATES ===")
    print(A.pivot(index="species",columns="scenario",values="n_supported").reindex(ORDER).to_string())
    print("\n=== Colobus original failures (loci) ===")
    print(failed[["group","locus_id","n_eligible_reads"]].to_string(index=False))

    for group in ("ALL","BLAST_QC_PASS","BLAST_QC_FAIL"):
        sub=A if group=="ALL" else S[S.group.eq(group)]
        if sub.empty: continue
        scenarios_plot=["Original","Flank_50","Flank_30","Locus_80","Locus_70",
                        "Reads_1","Moderate","Exploratory"]
        piv=sub.pivot(index="species",columns="scenario",values="rate_pct").reindex(ORDER)
        fig,ax=plt.subplots(figsize=(13,5))
        x=np.arange(len(piv))
        n=len(scenarios_plot)
        width=.85/n
        for i,sc in enumerate(scenarios_plot):
            ax.bar(x+(i-(n-1)/2)*width,piv.get(sc,pd.Series(index=piv.index,dtype=float)).to_numpy(),
                   width,label=sc)
        ax.set_xticks(x,[s.replace("_"," ") for s in piv.index],rotation=35,ha="right")
        ax.set_ylim(0,105)
        ax.set_ylabel("Nuclear-supported NUMTs (%)")
        ax.set_title(f"Nuclear threshold sensitivity: {group}")
        ax.legend(ncol=4,loc="upper center",bbox_to_anchor=(.5,-.27),frameon=False)
        fig.tight_layout()
        fig.savefig(out/f"threshold_sensitivity_{group}.png",dpi=220,bbox_inches="tight")
        plt.close(fig)
    print("\nSaved to",out)
    print("NOTE: identity and MAPQ stay fixed; no read-level chrM BLAST is rerun.")
    print("NOTE: lowering flank threshold to 30/50 bp in this script assesses those")
    print("      aligned bases within original 100-bp windows, not shorter new windows.")


if __name__=="__main__":
    main()
