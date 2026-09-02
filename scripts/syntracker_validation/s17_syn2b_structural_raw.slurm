#!/bin/bash
#SBATCH --job-name=syntracker_syn2b_raw
#SBATCH --partition=amd
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/data/syntracker_validation/logs/syn2b_structural_raw_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/data/syntracker_validation/logs/syn2b_structural_raw_%j.err
# Task 7 from HPC_BRIEFING.md: Syn2b structural channel on SynTracker cohorts, raw assemblies.
set -euo pipefail
BASE=/lustre1/g/aos_shihuang/data/syntracker_validation
ROOT=/lustre1/g/aos_shihuang/Syn2bANI-paper
PY=/group/aos_shihuang/conda/bin/python3

mkdir -p "$BASE/logs"

export RAYON_NUM_THREADS=1

$PY "$ROOT/scripts/syntracker_validation/08_syn2b_structural.py" \
    --assembly-dir "$BASE/assemblies" \
    --samples-dir  "$BASE/samples" \
    --out-dir      "$BASE/syn2b_structural_raw" \
    --syn2b        /lustre1/g/aos_shihuang/Syn2b/target/release/syn2b \
    --label        raw \
    --workers      16
