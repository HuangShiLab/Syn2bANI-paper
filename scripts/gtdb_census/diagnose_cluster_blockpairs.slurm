#!/bin/bash
#SBATCH --job-name=gtdb-blockqc
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/data/syn2b_census/logs/diagnose_cluster_blockpairs_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/data/syn2b_census/logs/diagnose_cluster_blockpairs_%j.err
set -euo pipefail
L=/lustre1/g/aos_shihuang
W=$L/data/syn2b_census
R=$L/Syn2b-paper/results/census/gtdb
cd "$L/Syn2b-paper"
python3 scripts/gtdb_census/diagnose_cluster_blockpairs.py \
  --table "$R/gtdb_r207_effective_normalized_dedup.tsv.gz" \
  --workdir "$W" \
  --cluster "Escherichia coli" \
  --block-size 500 \
  --output "$R/ecoli_blockpair_coverage.tsv" \
  --report "$R/ecoli_blockpair_coverage.json"
