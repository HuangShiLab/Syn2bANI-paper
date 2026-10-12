#!/bin/bash
#SBATCH --job-name=gtdb-qual-repair
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL
#SBATCH --output=/lustre1/g/aos_shihuang/data/syn2b_census/logs/gtdb_quality_sensitivity_repaired_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/data/syn2b_census/logs/gtdb_quality_sensitivity_repaired_%j.err
set -euo pipefail
L=/lustre1/g/aos_shihuang
cd "$L/Syn2b-paper"
mkdir -p "$L/Syn2b-paper/results/census/gtdb/effective_repaired"
python3 scripts/gtdb_census/gtdb_quality_sensitivity.py \
  --table "$L/Syn2b-paper/results/census/gtdb/gtdb_r207_effective_normalized_dedup.repaired.tsv.gz" \
  --quality "$L/data/syn2b_census/gtdb_genome_quality.tsv" \
  --output "$L/Syn2b-paper/results/census/gtdb/effective_repaired/effective_quality_sensitivity.tsv.tmp"
mv "$L/Syn2b-paper/results/census/gtdb/effective_repaired/effective_quality_sensitivity.tsv.tmp" \
   "$L/Syn2b-paper/results/census/gtdb/effective_repaired/effective_quality_sensitivity.tsv"
