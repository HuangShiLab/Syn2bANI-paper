#!/bin/bash
#SBATCH --job-name=ha_skani_final
#SBATCH --partition=intel
#SBATCH --array=0-44
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/ha_skani_final_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/ha_skani_final_%A_%a.err
# s8: skani dist for final high-ANI pairs.
set -euo pipefail
WORK=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k
export SYN2BANI_WORK=$WORK NSLICES=45
mkdir -p "$WORK/logs"
NEED=${NSLICES}
if [ "$SLURM_ARRAY_TASK_ID" -ge "$NEED" ]; then
    exit 0
fi
"$WORK/scripts/run_skani_high_ani_final_slice.sh" "$SLURM_ARRAY_TASK_ID"
