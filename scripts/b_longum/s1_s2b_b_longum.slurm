#!/bin/bash
#SBATCH --job-name=blongum_s2b
#SBATCH --partition=intel
#SBATCH --array=0-19
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=4:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/b_longum_abfA/logs/s2b_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/b_longum_abfA/logs/s2b_%A_%a.err
# s1: Syn2bANI all-vs-all for B. longum abfA case study (~180 genomes, ~16k pairs).
set -euo pipefail
WORK=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/b_longum_abfA
export B_LONGUM_WORK=$WORK NSLICES=20
mkdir -p "$WORK/logs"
# run_s2b_b_longum_slice.sh partitions the pair file into NSLICES equal
# chunks, so every task id below NSLICES has work; tasks above it exit.
if [ "$SLURM_ARRAY_TASK_ID" -ge "$NSLICES" ]; then
    exit 0
fi
"$WORK/scripts/run_s2b_b_longum_slice.sh" "$SLURM_ARRAY_TASK_ID"
