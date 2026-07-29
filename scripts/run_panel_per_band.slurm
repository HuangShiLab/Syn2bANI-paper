#!/bin/bash
#SBATCH --job-name=s2b_panel
#SBATCH --partition=amd
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=panel_per_band_%j.out
#SBATCH --error=panel_per_band_%j.err

set -euo pipefail
cd /lustre1/g/aos_shihuang/Syn2bANI-paper

python3 scripts/run_panel_per_band.py \
    --sample results/sample_anim_truth.tsv \
    --anim-dir anim_results \
    --strata results/gtdb_r207_100k_strata.tsv \
    --syn2bani /lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani \
    --outdir results/panel_by_band
