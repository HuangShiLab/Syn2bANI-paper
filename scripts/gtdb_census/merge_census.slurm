#!/bin/bash
#SBATCH --job-name=census-merge
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --mail-type=FAIL
#SBATCH --output=/lustre1/g/aos_shihuang/data/syn2b_census/logs/merge_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/data/syn2b_census/logs/merge_%j.err
set -euo pipefail
L=/lustre1/g/aos_shihuang
cd "$L/Syn2b-paper"
mkdir -p results/census/gtdb
python3 scripts/gtdb_census/merge_census.py \
  --workdir "$L/data/syn2b_census" \
  --out results/census/gtdb \
  --output-name gtdb_r207_within_species_sv.tsv.gz \
  --summary-sample-size 1000
echo GTDB-STRUCTURAL-MERGE-DONE
