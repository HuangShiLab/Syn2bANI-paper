#!/bin/bash
#SBATCH --job-name=g50k_fmh2000_invfrac
#SBATCH --partition=amd
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/syn2b_fmh2000_invfrac_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/syn2b_fmh2000_invfrac_%j.err
# FracMinHash scale 2000 (BcgI-sparse density) inverted-fraction pass over held-out 50k.
set -euo pipefail
ROOT=/lustre1/g/aos_shihuang/Syn2bANI-paper
WORK=$ROOT/results/gtdb50k
PY=/lustre1/g/aos_shihuang/tools/syn2b_py/bin/python

mkdir -p "$WORK/logs"

export RAYON_NUM_THREADS=1

$PY "$ROOT/scripts/run_syn2b_inverted_fraction.py" \
    --pairs      "$WORK/pairs_50k.tsv" \
    --genome-dir /lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all \
    --syn2b      /lustre1/g/aos_shihuang/Syn2b/target/release/syn2b \
    --mode       fracminhash --kmer 31 --scale 2000 \
    --tgt-dir    "$WORK/syn2b_tgts_cache_fmh2000" \
    --out        "$WORK/syn2b_inverted_fraction_50k_fmh2000.tsv" \
    --workers    16
