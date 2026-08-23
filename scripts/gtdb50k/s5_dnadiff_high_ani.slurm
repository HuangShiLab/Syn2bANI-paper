#!/bin/bash
#SBATCH --job-name=ha_dnadiff
#SBATCH --partition=intel
#SBATCH --array=0-79
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/ha_dnadiff_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/ha_dnadiff_%A_%a.err
# s5: dnadiff ANIm truth for high-ANI candidate pairs (8,000 pairs, 80 slices).
set -euo pipefail
WORK=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k
export SYN2BANI_WORK=$WORK NSLICES=80
mkdir -p "$WORK/logs"
NPAIRS=$(($(wc -l < "$WORK/high_ani_candidates.tsv") - 1))
NEED=$(( (NPAIRS + 99) / 100 ))
if [ "$SLURM_ARRAY_TASK_ID" -ge "$NEED" ]; then
    exit 0
fi
"$WORK/scripts/run_dnadiff_high_ani_slice.sh" "$SLURM_ARRAY_TASK_ID"
