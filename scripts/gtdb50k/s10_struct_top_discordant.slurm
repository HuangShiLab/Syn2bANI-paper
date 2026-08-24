#!/bin/bash
#SBATCH --job-name=struct_top_disc
#SBATCH --partition=intel
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/struct_top_discordant_%A.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k/logs/struct_top_discordant_%A.err
# s10: syn2bani struct for top discordant high-ANI GTDB pairs.
set -euo pipefail

WORK=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k
export SYN2BANI_WORK=$WORK
mkdir -p "$WORK/logs"

"$WORK/scripts/run_struct_top_discordant.sh" "$WORK/top_discordant_high_ani_pairs.tsv"
