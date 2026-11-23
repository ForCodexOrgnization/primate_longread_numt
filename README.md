# Primate long-read NUMT validation — Bouchet pilot

Standalone, **non-destructive** pilot; no changes to existing variant calling or filtering.

## Provenance
- Bouchet Slurm partition `day`, genome paths, manager paths: `primate_run_manager/config/bouchet.sh`.
- NUMT besthit path, conda env and tool modules: `primate_mt_variant_calling/nextflow.config`.
- Candidate besthit columns: `primate_mt_variant_calling/numt_detection/build_numt_besthit_table.sh`.
- The metadata provides **sample-level** platform mixtures; `run_accessions` are not labeled individually. ENA run-level verification is mandatory.

## Setup on Bouchet

```bash
# Upload/extract this directory to a project/scratch location; do not run large jobs in $HOME.
cd /path/to/primate_longread_numt
python3 scripts/prepare_metadata.py --metadata data/primate_metadata_master_FINAL_with_missing_Zoonomia.tsv --outdir data
python3 scripts/resolve_ena_runs.py --candidate-runs data/candidate_runs.tsv --out data/resolved_runs.tsv
# Inspect data/resolved_runs.tsv: read_type must be PACBIO or ONT.
# Download only verified long-read FASTQ files (or use existing local long-read FASTQ files).
# Edit data/local_reads.tsv with one row per run: run_accession, local_path, read_preset.
source config/bouchet.env
python3 scripts/build_manifest.py --resolved data/resolved_runs.tsv \
 --reference-dir "$WHOLE_REF_DIR" --besthit-dir "$BESTHIT_DIR" \
 --local-reads data/local_reads.tsv --out data/ready_manifest.tsv
mkdir -p logs
N=$(($(wc -l < data/ready_manifest.tsv)-1)); test "$N" -gt 0
# Pilot one run first
sbatch --array=1 slurm/validate_array.sbatch
# Full: sbatch --array=1-${N}%3 slurm/validate_array.sbatch
```

## Environment requirements
Existing `numt_env` must supply `python` + `pysam`, `blastn`, `makeblastdb` and **minimap2**. Verify first; `minimap2` is not guaranteed by the old repo. `samtools` is loaded from Bouchet module. Do not install packages into an existing shared conda environment without coordination; clone it or create a new env if necessary.

## Output meaning / limitations
`nuclear_support.tsv`: independent long-read coverage/spanning support for reference-encoded candidate nuclear loci. `reference_numt.fa` and `reference_numt_to_chrM.tsv`: **reference sequence** comparison to mtDNA, **not** long-read consensus or independent validation of chrM boundaries.

**Critical:** Spanning the reference NUMT region is necessary supportive evidence but is not a complete validation of the mtDNA origin or correct NUMT sequence. Full orthogonal validation requires long-read consensus/local assembly of the insertion, checking nuclear flanks, and aligning that independent sequence to sample-matched chrM; that extended step is **not implemented** in this pilot. Long-read discordance should be interpreted with reference bias, coverage, structural variation and sample identity in mind. Do not interpret 0 spanning reads as a false positive without adequate spanning power. Multi-block / inverted / circular chrM alignments must be retained as multiple blocks rather than collapsing to one interval.

## Known path caveats
The repo paths are the actual configured strings, not verified active Bouchet filesystem contents. In particular `Ref_whole` file naming, location of `numt_vs_chrM.besthit.tsv`, and source long-read FASTQs must be confirmed on Bouchet; `build_manifest.py` refuses missing files. FASTQ accessions need ENA verified per-run; PacBio subreads vs HiFi should use different presets (`map-pb` vs `map-hifi`).

## Future extension
1. Group read support by **biological sample** and collapse overlapping runs to avoid treating runs as independent samples.
2. Extract/reassemble actual NUMT insertions from anchored long reads, refine junctions.
3. BLAST assembled NUMTs against sample-specific consensus chrM, preserve multi-hit alignments and calculate reciprocal overlap / breakpoint errors.
4. Independently discover missed NUMTs; benchmark only assessable events.

## Streaming pilot: aria2c FASTQ downloads

The streaming Slurm workflow (`scripts/submit_batch.py` -> `slurm/download_align_validate.sbatch`) uses `aria2c` for parallel HTTPS transfers from ENA. The download is staged at `data/longread_download_stage/<sample>/<run>/`, with a persistent `.aria2` sidecar for segmented resumption after an interrupted or time-limited job. Files are gzip-checked before promotion, and staged FASTQs are removed only after successful alignment and validation.

On Bouchet, verify `aria2c --version` in the job environment. Set `ARIA2_MODULE` or `ARIA2_BIN` in `config/bouchet.env` if needed. The pipeline deliberately fails with an actionable error if aria2 is unavailable. Use `ARIA2_CONNECTIONS=4` initially; larger values may not help on a shared or throttled ENA connection. `DOWNLOAD_MAX_SECONDS` defaults to 72000 (20 hours), preserving download progress ahead of the 23-hour partition walltime. Rerun the same manifest to resume.

**Migration:** legacy curl `.partial` files without an aria2 `.aria2` metadata sidecar are renamed to `.legacy_curl_partial.<timestamp>` rather than being mixed with aria2 segments. Previously failed gzip files stay quarantined and must not be used for alignment. Inspect and delete quarantined files after confirming they are no longer needed; these can occupy tens of GB.

## ONT-read-to-chrM homology check (no new alignment)

After `slurm/validate_chrm_candidates.sbatch` produces the nuclear support
audit, run:

```bash
mkdir -p logs
sbatch slurm/validate_read_chrm_homology.sbatch
```

Defaults: `SRS4115496` / `SRR26729452` / `Pongo_abelii`. Override
`SAMPLE`, `RUN`, `SPECIES` with `sbatch --export=ALL,SAMPLE=...,RUN=...,SPECIES=...`
for another completed sample.

The script reads the *already generated* `chrM_candidate_validation` locus/read
evidence, selects at most 40 reads per locus with strict nuclear support, extracts
the CIGAR-mapped read bases corresponding to the nuclear candidate (including
insertions), and BLASTs those read-derived sequences against chrM in the
species-specific whole-genome FASTA. It neither downloads FASTQs nor remaps BAMs.
The default exploratory chrM thresholds are alignment >=75 bp and identity
>=75%, with an expected-site window expanded by 250 bp. These are provisional
screening cutoffs, not a definitive NUMT pathogenicity or origin call.
An expected-site hit is a BLAST alignment overlapping the expanded original
chrM interval; it does not prove exact boundaries. At least 2 distinct reads
must pass for the locus-level support label.

Outputs under `results/<sample>/<run>/longread_chrM_homology/`:

- `nuclear_mapped_read_segments.fa`: selected read-derived NUMT segments.
- `read_segments_vs_chrM.blast.tsv`: all BLAST HSPs.
- `per_read_chrM_homology.tsv`: best qualifying BLAST hit and expected-site hit.
- `per_locus_chrM_homology.tsv`: read counts and group (BLAST_QC_PASS / FAIL).

The analysis is **conditional on nuclear support** and uses short-read-derived
chrM coordinates for expected-site checking. It is supplementary sequence
homology evidence, not de novo assembly or proof of insertion age. Read-derived
homology can also reflect reference and mapping biases. Inspect the BLAST HSPs,
especially when a locus has a very short original chrM hit.

## Representative cross-species NUMT validation pilot (16 samples, 16 species)

The checked-in `data/pilot_run_selection.tsv` represents **16 biological
samples / 16 species / 22 ENA runs**. Six samples have planned ONT/PacBio
cross-platform pairs; ten add taxonomic diversity. This is a *candidate*
selection, not a claim that all runs are executable yet: reference,
best-hit file, FASTQ access, and particularly PacBio CLR/HiFi type need
verification. **Do not assign PacBio presets from instrument name alone.**

### Prepare on Bouchet

```bash
cd /nfs/roberts/scratch/pi_njl27/lt692/primate_longread_numt
python3 scripts/select_pilot_runs.py --resolved data/resolved_runs.tsv \
  --selection data/pilot_run_selection.tsv --out data/pilot_resolved_runs.tsv
python3 scripts/generate_presets.py --resolved data/pilot_resolved_runs.tsv \
  --out data/pilot_presets.tsv
# IMPORTANT: manually verify and replace all PacBio REVIEW with map-pb
# for CLR/subreads, or map-hifi for confirmed CCS/HiFi.
source config/bouchet.env
python3 scripts/build_remote_manifest.py \
  --resolved data/pilot_resolved_runs.tsv \
  --reference-dir "$WHOLE_REF_DIR" --besthit-dir "$BESTHIT_DIR" \
  --presets data/pilot_presets.tsv --out data/pilot_remote_manifest.tsv
cat data/pilot_remote_manifest.skipped.tsv
```

Only READY runs appear in `data/pilot_remote_manifest.tsv`.
Investigate skipped runs rather than asserting they completed. Missing
best-hit files may require generating them using the existing NUMT pipeline;
do not silently substitute data from another sample.

### Run alignments in bounded batches

```bash
python3 scripts/submit_batch.py --manifest data/pilot_remote_manifest.tsv \
  --unit run --batch-size 2 --concurrent 1 --dry-run
# Remove --dry-run only when ready. Do not launch duplicate runs while
# an earlier array is still active; rerun the command after job completion.
```

A completed run retains aligned BAM/BAI and `complete.flag`, even after
staged FASTQs are deleted. The submitter skips completed alignments.

### Validate completed alignments (no FASTQ transfer or realignment)

```bash
python3 scripts/submit_validation_batch.py \
  --manifest data/pilot_remote_manifest.tsv --batch-size 4 --concurrent 2 --dry-run
# Remove --dry-run to submit a Slurm array.
```

The validation submitter uses only completed alignment runs, and skips
runs with both existing validation outputs. Do not run multiple validation
submission batches concurrently for the same unfinished runs.

The new `slurm/validate_pilot_array.sbatch` runs, per sample/run:

1. `validate_chrm_candidates.py` (nuclear CIGAR and sequence support);
2. `validate_read_chrm_homology.py` (nuclear-mapped read fragments BLASTed
   against that species' chrM, including expected-locus homology).

For the first validated Pongo sample both existing results are reused.

### Summarize by sample, species and sequencing technology

```bash
python3 scripts/aggregate_representative_pilot.py \
  --manifest data/pilot_remote_manifest.tsv
```

Outputs in `results/representative_pilot_summary/`:

- `run_summary.tsv`: per-run status and counts;
- `run_locus_evidence.tsv`: per-locus per-run and platform evidence;
- `sample_species_summary.tsv`: non-duplicated sample/species counts,
  and validation status.

**Interpretation:** the numerator for both nuclear and chrM evidence is
`n_both_supported`; the denominator for the targeted analysis is the
number of chrM-hit loci, preferably stratified by `BLAST_QC_PASS` and
`BLAST_QC_FAIL`. Multiple runs from a biological sample must not be
treated as separate independent samples. The validation relies on
reference-guided nuclear coordinates and predicted chrM intervals; it
is targeted, supportive evidence, not de novo NUMT discovery.


## Expand from selected pilot to ALL available cohort long-read runs

The original `prepare_metadata.py` default produces *paired Illumina+long-read*
samples. It does **not** cover long-read-only samples. To include all cohort
samples with ONT or PacBio and keep the original pilot files untouched:

```bash
cd /nfs/roberts/scratch/pi_njl27/lt692/primate_longread_numt
python3 scripts/prepare_metadata.py \
  --metadata data/primate_metadata_master_FINAL_with_missing_Zoonomia.tsv \
  --outdir data --include-longread-only

# This can involve many ENA requests. Cached completed lookups are reused.
python3 scripts/resolve_ena_runs.py \
  --candidate-runs data/all_longread_candidate_runs.tsv \
  --out data/all_longread_resolved_runs.tsv

# Exclude the 22 explicitly selected pilot run accessions, but include any
# additional runs on pilot samples. Flag missing FASTQs and ENA errors.
python3 scripts/select_remaining_longread_runs.py \
  --resolved data/all_longread_resolved_runs.tsv \
  --pilot-selection data/pilot_run_selection.tsv \
  --out data/all_longread_nonpilot_resolved.tsv

# Review PacBio data type individually: map-pb (CLR/subreads) or
# map-hifi (CCS/HiFi). REVIEW is intentionally NOT allowed to run.
python3 scripts/generate_presets.py \
  --resolved data/all_longread_nonpilot_resolved.tsv \
  --out data/all_longread_nonpilot_presets.tsv
# Manually edit PacBio REVIEW entries in this TSV.

source config/bouchet.env
python3 scripts/build_remote_manifest.py \
  --resolved data/all_longread_nonpilot_resolved.tsv \
  --reference-dir "$WHOLE_REF_DIR" --besthit-dir "$BESTHIT_DIR" \
  --presets data/all_longread_nonpilot_presets.tsv \
  --out data/all_longread_nonpilot_manifest.tsv

# Inspect missing references, missing per-sample NUMT besthit, missing
# verified FASTQs, and unreviewed PacBio presets BEFORE submission.
cat data/all_longread_nonpilot_manifest.skipped.tsv

mkdir -p logs
python3 scripts/submit_batch.py \
  --manifest data/all_longread_nonpilot_manifest.tsv \
  --batch-size 2 --concurrent 1 --dry-run
```

Remove `--dry-run` only after confirming the manifest and that no corresponding
Slurm jobs are active. `submit_batch.py` skips completed alignment runs but
**does not track pending/running jobs**; running repeated submit commands without
checking `squeue -u "$USER"` can produce duplicate attempts. Proceed in bounded
batches, then validate completed alignments:

```bash
python3 scripts/submit_validation_batch.py \
  --manifest data/all_longread_nonpilot_manifest.tsv \
  --batch-size 4 --concurrent 2 --dry-run

# Remove --dry-run to submit only after reviewing the pending jobs.
# Aggregate new cohort separately, preserving existing pilot summaries:
python3 scripts/aggregate_representative_pilot.py \
  --manifest data/all_longread_nonpilot_manifest.tsv \
  --outdir results/all_longread_nonpilot_summary
```

Note that `build_remote_manifest.py` checks for a single species-named reference
and sample-matched NUMT besthit. Any missing input is skipped, **not** silently
replaced with a different species' reference. A skipped run is not a validated
run; investigate `*.skipped.tsv`. ENA run-level platform lookup is essential:
the original metadata lists platforms at the *sample* level, not the run level.
