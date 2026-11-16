#!/usr/bin/env bash
#SBATCH --job-name=primate_numt_validation
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=06:00:00
#SBATCH --output=logs/primate_validation_%A_%a.log
set -Eeuo pipefail
ROOT="${LR_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-}}"
[[ -n "$ROOT" && -s "$ROOT/config/bouchet.env" ]] || { echo "Invalid project root: $ROOT" >&2; exit 10; }
ROOT="$(cd "$ROOT" && pwd)"
source "$ROOT/config/bouchet.env"
: "${VALIDATION_MANIFEST:?Missing VALIDATION_MANIFEST}"
: "${SLURM_ARRAY_TASK_ID:?Not an array task}"
module load miniconda/24.11.3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$NUMT_ENV"
command -v blastn >/dev/null || { echo "blastn unavailable" >&2; exit 11; }
command -v makeblastdb >/dev/null || { echo "makeblastdb unavailable" >&2; exit 11; }
python -c 'import pysam'
python "$ROOT/scripts/validate_pilot_run.py" \
  --root "$ROOT" --manifest "$VALIDATION_MANIFEST" \
  --task-id "$SLURM_ARRAY_TASK_ID"
