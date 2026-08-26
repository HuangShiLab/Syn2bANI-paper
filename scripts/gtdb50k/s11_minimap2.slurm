#!/bin/bash
#SBATCH --job-name=g50k_mm2
#SBATCH --partition=intel
#SBATCH --array=0-189
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=2:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/minimap2_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/minimap2_%A_%a.err
# s11: minimap2 PAF for the 50k GTDB held-out benchmark. 190 slices of
#      ~230 pairs each; idempotent, safe to resubmit.
set -euo pipefail
WORK=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k
# NSLICES is overridable at submit time (sbatch --export=ALL,NSLICES=50
# --array=0-49) to fit account submit-slot limits.
NSLICES=${NSLICES:-190}
export GTDB50K_WORK=$WORK NSLICES
mkdir -p "$WORK/logs"
# The slice script partitions pairs into NSLICES equal chunks; every task
# id below NSLICES has work.
if [ "$SLURM_ARRAY_TASK_ID" -ge "$NSLICES" ]; then
    exit 0
fi
"$WORK/scripts/run_minimap2_slice.sh" "$SLURM_ARRAY_TASK_ID"
