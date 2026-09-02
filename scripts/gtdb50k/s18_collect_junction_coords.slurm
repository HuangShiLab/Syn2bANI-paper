#!/bin/bash
#SBATCH --job-name=g50k_junction_coords
#SBATCH --partition=amd
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/junction_coords_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/junction_coords_%j.err
# Task 9 from HPC_BRIEFING.md: collect Syn2b and dnadiff junction coordinates.
set -euo pipefail
ROOT=/lustre1/g/aos_shihuang/Syn2bANI-paper
WORK=$ROOT/results/gtdb50k
PY=/group/aos_shihuang/conda/bin/python3

mkdir -p "$WORK/logs"

$PY "$ROOT/scripts/gtdb50k/collect_junction_coordinates.py" \
    --pairs     "$WORK/pairs_50k.tsv" \
    --tgt-cache "$WORK/syn2b_tgts_cache" \
    --dnadiff   "$WORK/out" \
    --outdir    "$WORK/junction_coords" \
    --workers   32
