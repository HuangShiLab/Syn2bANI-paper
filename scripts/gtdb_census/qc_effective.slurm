#!/bin/bash
#SBATCH --job-name=gtdb-qc
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/data/syn2b_census/logs/gtdb_effective_qc_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/data/syn2b_census/logs/gtdb_effective_qc_%j.err
set -euo pipefail
L=/lustre1/g/aos_shihuang
cd "$L/Syn2b-paper"
python3 scripts/gtdb_census/qc_normalized_census.py \
  --table "$L/Syn2b-paper/results/census/gtdb/gtdb_r207_effective_normalized_dedup.tsv.gz" \
  --workdir "$L/data/syn2b_census" \
  --out "$L/Syn2b-paper/results/census/gtdb"
echo GTDB-EFFECTIVE-QC-DONE
