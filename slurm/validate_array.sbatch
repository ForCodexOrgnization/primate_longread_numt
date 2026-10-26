#!/usr/bin/env bash
#SBATCH --job-name=lr_numt_val
#SBATCH --partition=day
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=23:00:00
#SBATCH --output=logs/lr_numt_%A_%a.log
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$ROOT/config/bouchet.env"
module purge
module load SAMtools/1.21-GCC-13.3.0
module load miniconda/24.11.3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$NUMT_ENV"
command -v minimap2 >/dev/null || { echo 'minimap2 missing: install to NUMT_ENV' >&2; exit 1; }
python -c 'import pysam' || { echo 'pysam missing' >&2; exit 1; }
ROW=$(awk -F '\t' -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1{print;exit}' "$ROOT/data/ready_manifest.tsv")
[[ -n "$ROW" ]] || { echo 'no manifest row' >&2; exit 1; }
IFS=$'\t' read -r SAMPLE SPECIES RUN PLATFORM PRESET READS REF BESTHIT <<< "$ROW"
OUT="$ROOT/results/${SAMPLE}/${RUN}"
mkdir -p "$OUT"
[[ -s "${REF}.fai" ]] || { echo "Missing reference index: ${REF}.fai" >&2; exit 1; }
# Keep local output on project storage; avoid destructive cleanup.
if [[ ! -s "$OUT/aligned.bam.bai" ]]; then
 minimap2 -t 12 -ax "$PRESET" "$REF" "$READS" | samtools sort -@ 4 -m 2G -o "$OUT/aligned.bam" -
 samtools index "$OUT/aligned.bam"
fi
python "$ROOT/scripts/validate_loci.py" --bam "$OUT/aligned.bam" --reference "$REF" --besthit "$BESTHIT" --sample "$SAMPLE" --run "$RUN" --outdir "$OUT"
# Compare reference NUMT segments to species mtDNA with BLAST; this is a REFERENCE comparison,
# not independent long-read sequence reconstruction.
samtools faidx "$REF" chrM > "$OUT/chrM.fa"
makeblastdb -in "$OUT/chrM.fa" -dbtype nucl >/dev/null
blastn -query "$OUT/reference_numt.fa" -db "$OUT/chrM.fa" -outfmt '6 qseqid sseqid pident length qstart qend sstart send bitscore' > "$OUT/reference_numt_to_chrM.tsv"
touch "$OUT/complete.flag"
