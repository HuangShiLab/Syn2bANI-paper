#!/bin/bash
#SBATCH --job-name=mv_collect
#SBATCH --partition=intel
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/collect_%A.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/collect_%A.err
# s10: merge everything into $WORK/collect/*.tsv
set -euo pipefail
SCRIPTS_HOME=${SCRIPTS_HOME:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/scripts}
source "$SCRIPTS_HOME/common.sh"
mkdir -p "$WORK/logs"
$PY "$SCRIPTS/collect.py" "$WORK"
echo "[collect] outputs:"
ls -la "$WORK/collect/"
