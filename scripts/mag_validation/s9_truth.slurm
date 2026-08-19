#!/bin/bash
#SBATCH --job-name=mv_dnadiff
#SBATCH --partition=intel
#SBATCH --array=0-199
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/truth_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/truth_%A_%a.err
# s9: dnadiff ANIm truth, MAG vs majority source genome. Chunked array:
#     each task handles DNADIFF_CHUNK pairs; out-of-range tasks exit instantly.
set -euo pipefail
SCRIPTS_HOME=${SCRIPTS_HOME:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/scripts}
source "$SCRIPTS_HOME/common.sh"
mkdir -p "$WORK/logs"
NPAIRS=$(($(wc -l < "$WORK/pairs/pairs_anchor.tsv") - 1))
NEED=$(( (NPAIRS + DNADIFF_CHUNK - 1) / DNADIFF_CHUNK ))
if [ "$SLURM_ARRAY_TASK_ID" -ge "$NEED" ]; then
    exit 0
fi
"$SCRIPTS/run_dnadiff_chunk.sh" "$SLURM_ARRAY_TASK_ID"
