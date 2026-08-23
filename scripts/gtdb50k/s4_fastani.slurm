#!/bin/bash
#SBATCH --job-name=g50k_fastani
#SBATCH --partition=intel
#SBATCH --array=0-189
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/fastani_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/fastani_%A_%a.err
# s4: FastANI on the 50k GTDB held-out pairs (190 slices of ~230 pairs).
#     Idempotent per-slice; submit in waves of <=25 (submit quota).
set -euo pipefail
WORK=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k
export GTDB50K_WORK=$WORK NSLICES=190
mkdir -p "$WORK/logs"
NPAIRS=$(($(wc -l < "$WORK/pairs_50k.tsv") - 1))
NEED=$(( (NPAIRS + 228) / 229 ))
if [ "$SLURM_ARRAY_TASK_ID" -ge "$NEED" ]; then
    exit 0
fi
ROWS_OUT="$WORK/fastani_rows/slice_${SLURM_ARRAY_TASK_ID}.tsv"
if [ -s "$ROWS_OUT" ]; then
    exit 0
fi
"$WORK/scripts/run_fastani_slice.sh" "$SLURM_ARRAY_TASK_ID"
