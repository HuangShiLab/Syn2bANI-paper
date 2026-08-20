#!/bin/bash
#SBATCH --job-name=mv_checkm2
#SBATCH --partition=intel
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/checkm2_%A.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/checkm2_%A.err
# s4: stage all bins >=100kbp into bins_all/ (headers rewritten), compute bin
#     stats, then CheckM2 predict on the whole cohort.
set -euo pipefail
SCRIPTS_HOME=${SCRIPTS_HOME:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/scripts}
source "$SCRIPTS_HOME/common.sh"
mkdir -p "$WORK/logs" "$WORK/bins_all"

$PY "$SCRIPTS/stage_bins.py" "$WORK/bins" "$WORK/bins_all" "$WORK/bins_all/bin_stats.tsv" "$MIN_BIN_BP"

if [ -s "$WORK/checkm2/quality_report.tsv" ]; then
    echo "[checkm2] quality_report.tsv exists, skipping"
    exit 0
fi
export PYTHONNOUSERSITE=1                  # ~/.local numpy shadows the env and breaks tensorflow
export PATH="$(dirname "$CHECKM2"):$PATH"   # diamond lives in the same env
$CHECKM2 predict \
    --input "$WORK/bins_all" -x fa \
    --output-directory "$WORK/checkm2" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --database_path "$CHECKM2DB" --force
echo "[checkm2] done"
