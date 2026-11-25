#!/usr/bin/env bash
#SBATCH --job-name=lr_numt_manager
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=logs/cohort_manager_%j.log
set -Eeuo pipefail
ROOT="${LR_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-}}"
: "${ROOT:?Missing LR_PROJECT_ROOT}"
: "${COHORT_MANIFEST:?Missing COHORT_MANIFEST}"
cd "$ROOT"
ARGS=(--manifest "$COHORT_MANIFEST" --submit
  --min-free-gib "${COHORT_MIN_FREE_GIB:-500}"
  --max-align "${COHORT_MAX_ALIGN:-1}"
  --max-validate "${COHORT_MAX_VALIDATE:-1}"
  --max-merge "${COHORT_MAX_MERGE:-1}"
  --max-attempts "${COHORT_MAX_ATTEMPTS:-5}")
if [[ "${COHORT_CLEANUP:-0}" == 1 ]]; then ARGS+=(--cleanup-if-validated); fi
if [[ -e data/cohort_manager/STOP ]]; then echo "STOP marker present; not scheduling."; exit 0; fi
python3 scripts/cohort_manager.py "${ARGS[@]}"
# Another controller can start only when this tick has finished successfully.
# Stop when every run and merge stage has completed.
if [[ -e data/cohort_manager/STOP ]]; then exit 0; fi
if python3 - <<'PY'
import json,sys
from pathlib import Path
p=Path("data/cohort_manager/state.json")
if not p.is_file():sys.exit(1)
j=json.loads(p.read_text()).get("jobs",{})
sys.exit(0 if j and all(v.get("status")=="DONE" for v in j.values()) else 1)
PY
then echo "ALL STAGES COMPLETE"; exit 0; fi
NUMBER="${COHORT_TICK_NO:-1}"
MAX="${COHORT_MAX_TICKS:-10000}"
if (( NUMBER >= MAX )); then echo "Maximum manager ticks reached; inspect status"; exit 0; fi
NEXT=$((NUMBER+1))
sbatch --parsable --dependency="afterok:${SLURM_JOB_ID:?}" --begin=now+5minutes \
  --export="ALL,LR_PROJECT_ROOT=$ROOT,COHORT_MANIFEST=$COHORT_MANIFEST,COHORT_MIN_FREE_GIB=${COHORT_MIN_FREE_GIB:-500},COHORT_MAX_ALIGN=${COHORT_MAX_ALIGN:-1},COHORT_MAX_VALIDATE=${COHORT_MAX_VALIDATE:-1},COHORT_MAX_MERGE=${COHORT_MAX_MERGE:-1},COHORT_MAX_ATTEMPTS=${COHORT_MAX_ATTEMPTS:-5},COHORT_CLEANUP=${COHORT_CLEANUP:-0},COHORT_TICK_NO=$NEXT,COHORT_MAX_TICKS=$MAX" \
  "$ROOT/slurm/cohort_manager_tick.sbatch"
echo "Next manager tick submitted: $NEXT"
