#!/usr/bin/env python3
"""Build a non-destructive non-pilot long-read run list from ENA-verified metadata.

The output contains all verified ONT/PacBio runs not explicitly selected in
pilot_run_selection.tsv. Runs on previously selected biological samples are
included if their run accession was not selected for the pilot.

No run is scheduled here. Missing ENA FASTQs, unknown platforms and duplicate
run/sample conflicts are reported separately for review.
"""
import argparse
import csv
from collections import Counter
from pathlib import Path


def read(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        wr = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        wr.writeheader()
        wr.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--resolved", default="data/all_longread_resolved_runs.tsv")
    ap.add_argument("--pilot-selection", default="data/pilot_run_selection.tsv")
    ap.add_argument("--out", default="data/all_longread_nonpilot_resolved.tsv")
    args = ap.parse_args()
    rows = read(args.resolved)
    pilot = read(args.pilot_selection)
    pilot_runs = {r["run_accession"] for r in pilot}
    selected = []
    skipped = []
    seen = {}
    for r in rows:
        run = r.get("run_accession", "").strip()
        identity = (r.get("sample", ""), r.get("species", ""))
        if not run:
            skipped.append(dict(**r, reason="missing_run_accession"))
        elif run in seen:
            reason = "duplicate_run_conflicting_identity" if seen[run] != identity else "duplicate_run"
            skipped.append(dict(**r, reason=reason))
        else:
            seen[run] = identity
            reason = ("selected_pilot_run" if run in pilot_runs else
                      "non_longread_platform" if r.get("read_type") not in ("ONT", "PACBIO") else
                      "ena_lookup_error" if r.get("lookup_error", "").strip() else
                      "fastq_missing" if not r.get("fastq_ftp", "").strip() else "")
            if reason:
                skipped.append(dict(**r, reason=reason))
            else:
                selected.append(r)
    if not rows:
        raise SystemExit("Input resolved runs is empty")
    fields = list(rows[0])
    write(args.out, selected, fields)
    write(Path(args.out).with_suffix(".skipped.tsv"), skipped, fields + ["reason"])
    def ns(rs):
        return len({x.get("sample") for x in rs})
    print("Verified ENA candidates:",len(rows))
    print("Pilot-selected runs:",len(pilot_runs))
    print("Non-pilot runnable long-read runs:",len(selected),
          "samples:",ns(selected),"species:",len({x["species"] for x in selected}))
    print("Skipped:",dict(Counter(r["reason"] for r in skipped)))
    print("Output:",args.out)
    print("IMPORTANT: PacBio read presets still require per-run verification.")


if __name__ == "__main__":
    main()
