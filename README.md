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
