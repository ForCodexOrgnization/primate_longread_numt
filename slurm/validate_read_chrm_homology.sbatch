#!/usr/bin/env bash
#SBATCH --job-name=lr_numt_chrm_homology
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=06:00:00
#SBATCH --output=logs/lr_chrm_homology_%j.log
set -Eeuo pipefail

ROOT="${LR_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-}}"
[[ -n "$ROOT" && -s "$ROOT/config/bouchet.env" ]] || { echo "Invalid LR_PROJECT_ROOT: $ROOT" >&2; exit 10; }
ROOT="$(cd "$ROOT" && pwd)"
source "$ROOT/config/bouchet.env"

SAMPLE="${SAMPLE:-SRS4115496}"
RUN="${RUN:-SRR26729452}"
SPECIES="${SPECIES:-Pongo_abelii}"
BAM="${BAM:-$ROOT/results/$SAMPLE/$RUN/aligned.bam}"
REFERENCE="${REFERENCE:-$WHOLE_REF_DIR/$SPECIES.fa}"
NUCLEAR_VALIDATION="${NUCLEAR_VALIDATION:-$ROOT/results/$SAMPLE/$RUN/chrm_candidate_validation}"
OUTDIR="${OUTDIR:-$ROOT/results/$SAMPLE/$RUN/longread_chrM_homology}"

for path in "$BAM" "$BAM.bai" "$REFERENCE" "$REFERENCE.fai" \
 "$NUCLEAR_VALIDATION/chrM_candidate_locus_validation.tsv" \
 "$NUCLEAR_VALIDATION/chrM_candidate_read_evidence.tsv"; do
  [[ -s "$path" ]] || { echo "Missing input: $path" >&2; exit 23; }
done

module load miniconda/24.11.3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$NUMT_ENV"
command -v blastn >/dev/null || { echo "blastn not available in NUMT_ENV=$NUMT_ENV" >&2; exit 11; }
python -c 'import pysam'
python "$ROOT/scripts/validate_read_chrm_homology.py" \
  --bam "$BAM" --reference "$REFERENCE" \
  --nuclear-validation "$NUCLEAR_VALIDATION" \
  --outdir "$OUTDIR" --threads "${SLURM_CPUS_PER_TASK:-2}"
echo "[COMPLETE] chrM homology audit: $OUTDIR"
