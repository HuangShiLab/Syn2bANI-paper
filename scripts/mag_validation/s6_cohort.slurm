#!/bin/bash
#SBATCH --job-name=mv_cohort
#SBATCH --partition=intel
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/cohort_%A.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/logs/cohort_%A.err
# s6: merge bin stats + CheckM2 + assignment -> collect/bins.tsv, pairs manifests.
set -euo pipefail
SCRIPTS_HOME=${SCRIPTS_HOME:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/mag_validation/scripts}
source "$SCRIPTS_HOME/common.sh"
mkdir -p "$WORK/logs"
$PY "$SCRIPTS/build_cohort.py" "$WORK"
