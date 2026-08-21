#!/bin/bash
#SBATCH --job-name=g50k_dnadiff
#SBATCH --partition=intel
#SBATCH --array=0-189
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/dnadiff_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/dnadiff_%A_%a.err
# s2: dnadiff ANIm truth for the 50k GTDB held-out benchmark. 190 slices of
#     ~230 pairs each (~40-70 min per slice); idempotent, safe to resubmit.
set -euo pipefail
WORK=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k
export GTDB50K_WORK=$WORK GTDB50K_SCRIPTS=$WORK/scripts NSLICES=190
mkdir -p "$WORK/logs"
NPAIRS=$(($(wc -l < "$WORK/pairs_50k.tsv") - 1))
NEED=$(( (NPAIRS + 229) / 230 ))
if [ "$SLURM_ARRAY_TASK_ID" -ge "$NEED" ]; then
    exit 0
fi
"$WORK/scripts/run_dnadiff_slice.sh" "$SLURM_ARRAY_TASK_ID"
