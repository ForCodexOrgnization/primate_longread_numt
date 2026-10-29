# Representative pilot (16 samples / 16 species / 22 runs)

This is an OPTIONAL restricted-run workflow. It does not alter the full manifest or the existing variant-calling pipeline.

## Selection

- `data/pilot_run_selection.tsv`: 22 verified ENA long-read run accessions, 16 biological samples and 16 species.
- Stage `cross_platform`: 12 runs from 6 samples, each with PacBio and ONT.
- Stage `diversity`: 10 runs, each from a different additional species.
- `scripts/select_pilot_runs.py` checks run accession, sample, species, long-read platform and available FASTQ.
- No automatic substitution of other runs if an input is missing.
- File availability, actual read quality and transfer size still need checking on Bouchet.

## Bouchet commands

```bash
cd /nfs/roberts/scratch/pi_njl27/lt692/primate_longread_numt
git fetch origin
git switch feature/representative-pilot-selection-20261029
mkdir -p logs

python3 scripts/select_pilot_runs.py \
  --resolved data/resolved_runs.tsv \
  --selection data/pilot_run_selection.tsv \
  --out data/pilot_resolved_runs.tsv

# For the 22 pilot runs only. PacBio values remain REVIEW until verified.
python3 scripts/generate_presets.py \
  --resolved data/pilot_resolved_runs.tsv \
  --out data/pilot_presets.tsv

# Edit data/pilot_presets.tsv:
# for verified CCS/HiFi -> map-hifi; for CLR/subreads -> map-pb
# REVIEW rows are deliberately skipped.

source config/bouchet.env
python3 scripts/build_remote_manifest.py \
  --resolved data/pilot_resolved_runs.tsv \
  --reference-dir "$WHOLE_REF_DIR" \
  --besthit-dir "$BESTHIT_DIR" \
  --presets data/pilot_presets.tsv \
  --out data/pilot_remote_manifest.tsv

# Inspect pilot_remote_manifest.skipped.tsv. Do not submit if zero ready.
wc -l data/pilot_remote_manifest.tsv
cat data/pilot_remote_manifest.skipped.tsv

python3 scripts/submit_batch.py \
  --manifest data/pilot_remote_manifest.tsv \
  --unit run --batch-size 2 --concurrent 1 --dry-run

# After review, remove --dry-run:
python3 scripts/submit_batch.py \
  --manifest data/pilot_remote_manifest.tsv \
  --unit run --batch-size 2 --concurrent 1

# After jobs complete; sample counts are not molecule-deduplicated:
python3 scripts/aggregate_sample_support.py \
  --manifest data/pilot_remote_manifest.tsv \
  --results-root results --outdir results/pilot_aggregation
```

**IMPORTANT**: `submit_batch.py` currently selects the earliest unfinished manifest rows, and does not reserve active rows. Wait until a prior batch finishes before submitting another. This is NOT an automatic queue-draining workflow.

**IMPORTANT**: The existing nuclear-reference best-hit locations are required. If these are missing or ambiguous, the manifest rejects those rows rather than silently choosing another reference. Check `BESTHIT_DIR`.

**IMPORTANT**: PacBio instrument names alone cannot establish HiFi versus raw subreads. Inspect the ENA filenames and study metadata, and confirm read type before replacing REVIEW. Do not blindly assign all Sequel II/Revio runs `map-hifi`.

**IMPORTANT**: The pipeline validates nuclear spanning support for reference-encoded predicted NUMTs; it does not yet independently reconstruct NUMT insertions or calculate long-read chrM interval accuracy.
