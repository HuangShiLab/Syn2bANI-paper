#!/bin/bash
#SBATCH --job-name=g50k_bcgI_invfrac
#SBATCH --partition=amd
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/syn2b_bcgI_invfrac_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/syn2b_bcgI_invfrac_%j.err
# Task 3 from HPC_TODO.md: single-enzyme BcgI inverted-fraction pass over held-out 50k.
set -euo pipefail
ROOT=/lustre1/g/aos_shihuang/Syn2bANI-paper
WORK=$ROOT/results/gtdb50k
PY=/group/aos_shihuang/conda/bin/python3

mkdir -p "$WORK/logs"

export RAYON_NUM_THREADS=1

$PY "$ROOT/scripts/run_syn2b_inverted_fraction.py" \
    --pairs      "$WORK/pairs_50k.tsv" \
    --genome-dir /lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all \
    --syn2b      /lustre1/g/aos_shihuang/Syn2b/target/release/syn2b \
    --enzymes    BcgI \
    --tgt-dir    "$WORK/syn2b_tgts_cache_bcgI" \
    --out        "$WORK/syn2b_inverted_fraction_50k_bcgI.tsv" \
    --workers    16
