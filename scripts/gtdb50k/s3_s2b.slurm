#!/bin/bash
#SBATCH --job-name=g50k_s2b
#SBATCH --partition=intel
#SBATCH --array=0-189
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/s2b_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/s2b_%A_%a.err
# s3: syn2bani ani --verbose --calibrate on all 43,334 held-out pairs.
#     Same 190-slice layout as s2; idempotent, safe to resubmit.
set -euo pipefail
WORK=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k
export GTDB50K_WORK=$WORK NSLICES=190
mkdir -p "$WORK/logs"
"$WORK/scripts/run_s2b_slice.sh" "$SLURM_ARRAY_TASK_ID"
