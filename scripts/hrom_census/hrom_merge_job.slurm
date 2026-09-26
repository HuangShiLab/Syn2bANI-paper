#!/bin/bash
# Run after census_worker --drain reports zero tasks.
#SBATCH --job-name=hrom-merge
#SBATCH --partition=amd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/databases/HROM/census/logs/merge_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/databases/HROM/census/logs/merge_%j.err

set -euo pipefail
REPO=/lustre1/g/aos_shihuang/Syn2b-paper
WORKDIR=/lustre1/g/aos_shihuang/databases/HROM/census
OUT=/lustre1/g/aos_shihuang/Syn2b-paper/results/census/hrom
mkdir -p "$WORKDIR/logs" "$OUT"
cd "$REPO"
python3 scripts/gtdb_census/merge_census.py \
    --workdir "$WORKDIR" \
    --out "$OUT" \
    --output-name hrom_within_species_sv.tsv.gz
