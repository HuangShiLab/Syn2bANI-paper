#!/bin/bash
#SBATCH --job-name=closed_invfrac
#SBATCH --partition=amd
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/data/closed_inversions/logs/closed_invfrac_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/data/closed_inversions/logs/closed_invfrac_%j.err
# Task 7b: all-vs-all Syn2b inverted_fraction on closed/high-quality genomes.
set -euo pipefail
ROOT=/lustre1/g/aos_shihuang/Syn2bANI-paper
WORK=/lustre1/g/aos_shihuang/data/closed_inversions
PY=/group/aos_shihuang/conda/bin/python3

mkdir -p "$WORK/results" "$WORK/logs"

export RAYON_NUM_THREADS=1

$PY "$ROOT/scripts/run_syn2b_inverted_fraction.py" \
    --pairs      "$WORK/pairs_all_vs_all.tsv" \
    --genome-dir "$WORK/genomes" \
    --syn2b      /lustre1/g/aos_shihuang/Syn2b/target/release/syn2b \
    --tgt-dir    "$WORK/syn2b_tgts_cache" \
    --out        "$WORK/results/syn2b_inverted_fraction_closed.tsv" \
    --workers    16
