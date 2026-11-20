#!/usr/bin/env bash
#SBATCH --job-name=lr_numt_candidate_audit
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=06:00:00
#SBATCH --output=logs/lr_candidate_audit_%j.log
set -Eeuo pipefail
ROOT="${LR_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-}}"
[[ -n "$ROOT" && -s "$ROOT/config/bouchet.env" ]] || { echo "Project root invalid: $ROOT" >&2; exit 10; }
ROOT="$(cd "$ROOT" && pwd)"
source "$ROOT/config/bouchet.env"
SAMPLE="${SAMPLE:-SRS4115496}"
RUN="${RUN:-SRR26729452}"
SPECIES="${SPECIES:-Pongo_abelii}"
BAM="${BAM:-$ROOT/results/$SAMPLE/$RUN/aligned.bam}"
REFERENCE="${REFERENCE:-$WHOLE_REF_DIR/$SPECIES.fa}"
BESTHIT="${BESTHIT:-$BESTHIT_DIR/$SAMPLE.numt_vs_chrM.besthit.tsv}"
OUTDIR="${OUTDIR:-$ROOT/results/$SAMPLE/$RUN/chrm_candidate_validation}"
for path in "$BAM" "$BAM.bai" "$REFERENCE" "$REFERENCE.fai" "$BESTHIT"; do
 [[ -s "$path" ]] || { echo "Missing input: $path" >&2; exit 23; }
done
module load miniconda/24.11.3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$NUMT_ENV"
python "$ROOT/scripts/validate_chrm_candidates.py" \
  --bam "$BAM" --reference "$REFERENCE" --besthit "$BESTHIT" --outdir "$OUTDIR"
echo "COMPLETE: $OUTDIR"
