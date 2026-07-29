#!/bin/bash
#SBATCH --job-name=s2b_stra_ref
#SBATCH --partition=amd
#SBATCH --array=1-25
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=strata_ref_%A_%a.out
#SBATCH --error=strata_ref_%A_%a.err

set -euo pipefail

cd /lustre1/g/aos_shihuang/Syn2bANI-paper

IDX="$SLURM_ARRAY_TASK_ID"
CHUNK_TSV=$(printf "results/strata_ref_chunks/pairs_%03d.tsv" "$IDX")
OUT_STRATA=$(printf "results/strata_ref_chunks/strata_%03d.tsv" "$IDX")
OUT_ANI=$(printf "results/strata_ref_chunks/ani_%03d.tsv" "$IDX")

echo "Starting chunk $IDX at $(date)"
echo "Chunk file: $CHUNK_TSV"

python3 scripts/run_strata_by_reference.py run "$CHUNK_TSV" \
    --syn2bani /lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani \
    --strata-out "$OUT_STRATA" \
    --ani-out "$OUT_ANI" \
    --threads-per-run 4 \
    --max-workers 8

echo "Finished chunk $IDX at $(date)"
