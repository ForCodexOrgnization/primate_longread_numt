#!/usr/bin/env bash
#SBATCH --job-name=numt_target_merge
#SBATCH --partition=day
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/numt_target_merge_%A_%a.log
set -Eeuo pipefail
ROOT="${LR_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-}}"
: "${ROOT:?Missing project root}"
: "${NUMT_MERGE_BATCH:?Missing merge batch TSV}"
: "${NUMT_MERGE_MANIFESTS:?Missing input manifests}"
: "${NUMT_MERGE_OUT_ROOT:?Missing out root}"
module load SAMtools/1.21-GCC-13.3.0
ROW="$(awk -F '\t' -v n="${SLURM_ARRAY_TASK_ID:?}" 'NR==n+1{print;exit}' "$NUMT_MERGE_BATCH" | tr -d '\015')"
[[ -n "$ROW" ]] || { echo "Missing batch row" >&2; exit 20; }
IFS=$'\t' read -r SAMPLE PLATFORM <<< "$ROW"
[[ "$SAMPLE" =~ ^[A-Za-z0-9._-]+$ ]] || exit 21
[[ "$PLATFORM" == ONT || "$PLATFORM" == PACBIO ]] || exit 22
IFS=',' read -r -a MANIFESTS <<< "$NUMT_MERGE_MANIFESTS"
COMMAND=(python3 "$ROOT/scripts/merge_sample_platform_numt_bam.py"
  --sample "$SAMPLE" --platform "$PLATFORM" --out-root "$NUMT_MERGE_OUT_ROOT"
  --results-root "$ROOT/results" --threads "${SLURM_CPUS_PER_TASK:-4}")
for manifest in "${MANIFESTS[@]}"; do COMMAND+=(--manifest "$manifest"); done
if [[ "${NUMT_MERGE_ALLOW_PARTIAL:-0}" == 1 ]]; then COMMAND+=(--allow-partial); fi
"${COMMAND[@]}"
