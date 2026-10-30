#!/usr/bin/env bash
#SBATCH --job-name=lr_numt_stream
#SBATCH --partition=day
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=23:00:00
#SBATCH --output=logs/lr_stream_%A_%a.log
set -Eeuo pipefail
# Slurm copies the script into /var/spool/slurmd; BASH_SOURCE is NOT the repo path.
ROOT="${LR_PROJECT_ROOT:-${SLURM_SUBMIT_DIR:-}}"
[[ -n "$ROOT" && -s "$ROOT/config/bouchet.env" && -s "$ROOT/scripts/validate_loci.py" ]] || {
  echo "ERROR: invalid project root: $ROOT (set LR_PROJECT_ROOT)" >&2
  exit 10
}
ROOT="$(cd "$ROOT" && pwd)"
source "$ROOT/config/bouchet.env"
: "${BATCH_MANIFEST:?Missing BATCH_MANIFEST}"
: "${LR_OUT_ROOT:?Missing LR_OUT_ROOT}"
module purge
module load SAMtools/1.21-GCC-13.3.0
module load miniconda/24.11.3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$NUMT_ENV"
for cmd in curl minimap2 samtools blastn makeblastdb python gzip; do command -v "$cmd" >/dev/null || { echo "Missing $cmd" >&2; exit 1; }; done
python -c 'import pysam'
ROW="$(awk -F '\t' -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1{print;exit}' "$BATCH_MANIFEST")"
[[ -n "$ROW" ]] || exit 20
IFS=$'\t' read -r SAMPLE SPECIES RUN PLATFORM PRESET URLS REF BESTHIT <<< "$ROW"
[[ "$SAMPLE" =~ ^[A-Za-z0-9._-]+$ && "$RUN" =~ ^[A-Za-z0-9._-]+$ ]] || exit 21
[[ "$PRESET" == map-pb || "$PRESET" == map-hifi || "$PRESET" == map-ont || "$PRESET" == lr:hq ]] || exit 22
[[ -s "$REF" && -s "${REF}.fai" && -s "$BESTHIT" ]] || exit 23
OUT="$LR_OUT_ROOT/$SAMPLE/$RUN"
mkdir -p "$OUT"
if [[ -f "$OUT/complete.flag" && -s "$OUT/aligned.bam" && -s "$OUT/aligned.bam.bai" && -s "$OUT/nuclear_support.tsv" && -s "$OUT/run_provenance.json" ]]; then
 samtools quickcheck "$OUT/aligned.bam"
 echo "Already completed: $RUN"; exit 0
fi
SCRATCH_BASE="${SLURM_TMPDIR:-${TMPDIR:-}}"
[[ -n "$SCRATCH_BASE" && -d "$SCRATCH_BASE" ]] || { echo 'SLURM_TMPDIR or TMPDIR required' >&2; exit 24; }
WORK="$(mktemp -d "$SCRATCH_BASE/lrnumt_${RUN}.XXXXXXXX")"
cleanup() { [[ ! -d "${WORK:-/nonexistent}" ]] || rm -rf -- "$WORK"; }
trap cleanup EXIT
IFS=';' read -r -a URL_ARRAY <<< "$URLS"
[[ "${#URL_ARRAY[@]}" -gt 0 ]] || exit 25
FILES=()
i=0
for url in "${URL_ARRAY[@]}"; do
 [[ "$url" == ftp.sra.ebi.ac.uk/* || "$url" == https://ftp.sra.ebi.ac.uk/* ]] || exit 26
 [[ "$url" == https://* ]] || url="https://$url"
 i=$((i+1)); ext=.fastq.gz
 [[ "$url" == *.fastq || "$url" == *.fq ]] && ext=.fastq
 file="$WORK/read_${i}${ext}"
 curl --fail --location --retry 4 --retry-delay 10 --connect-timeout 30 -o "$file.partial" "$url"
 [[ -s "$file.partial" ]] || exit 27
 mv "$file.partial" "$file"
 [[ "$file" == *.gz ]] && gzip -t "$file"
 FILES+=("$file")
done
reads_stream() {
 for file in "${FILES[@]}"; do
  if [[ "$file" == *.gz ]]; then gzip -dc "$file"; else cat "$file"; fi
 done
}
rm -f "$OUT/complete.flag"
reads_stream | minimap2 -t 12 -ax "$PRESET" "$REF" - | samtools sort -@ 4 -m 2G -o "$OUT/aligned.bam.partial" -
samtools quickcheck "$OUT/aligned.bam.partial"
mv "$OUT/aligned.bam.partial" "$OUT/aligned.bam"
samtools index "$OUT/aligned.bam"
samtools idxstats "$OUT/aligned.bam" >/dev/null
python "$ROOT/scripts/validate_loci.py" --bam "$OUT/aligned.bam" --reference "$REF" --besthit "$BESTHIT" --sample "$SAMPLE" --run "$RUN" --outdir "$OUT"
[[ -s "$OUT/nuclear_support.tsv" && -s "$OUT/reference_numt.fa" ]] || exit 28
samtools faidx "$REF" chrM > "$OUT/chrM.fa"
makeblastdb -in "$OUT/chrM.fa" -dbtype nucl >/dev/null
blastn -query "$OUT/reference_numt.fa" -db "$OUT/chrM.fa" -outfmt '6 qseqid sseqid pident length qstart qend sstart send bitscore' > "$OUT/reference_numt_to_chrM.tsv"
python - "$OUT" "$RUN" "$PRESET" "${#FILES[@]}" <<'PY'
import datetime,json,pathlib,sys
out,run,preset,n=sys.argv[1:]
(pathlib.Path(out)/'run_provenance.json').write_text(json.dumps(dict(run=run,preset=preset,downloaded_files=int(n),completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),indent=2))
PY
[[ -s "$OUT/aligned.bam.bai" && -s "$OUT/run_provenance.json" ]] || exit 29
touch "$OUT/complete.flag"
echo "COMPLETE: $RUN; staged FASTQs cleaned at job exit"
