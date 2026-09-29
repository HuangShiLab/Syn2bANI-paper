#!/bin/bash
#SBATCH --job-name=hrom-quality
#SBATCH --partition=amd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/databases/HROM/census/logs/quality_stratified_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/databases/HROM/census/logs/quality_stratified_%j.err

set -euo pipefail
REPO=/lustre1/g/aos_shihuang/Syn2bANI-paper
CENSUS=/lustre1/g/aos_shihuang/Syn2b-paper/results/census/hrom/hrom_within_species_sv_ani.tsv.gz
WORKDIR=/lustre1/g/aos_shihuang/databases/HROM/census
METADATA=/lustre1/g/aos_shihuang/databases/HROM/data/genome_catalog/HROM_representative_genome_metadata.tsv
GENOME_QUALITY=/lustre1/g/aos_shihuang/databases/HROM/data/genome_catalog/HROM_Conspecific-genomes-metadata.tsv
OUT=/lustre1/g/aos_shihuang/Syn2b-paper/results/census/hrom/quality_stratified
mkdir -p "$WORKDIR/logs" "$OUT"
cd "$REPO"
python3 scripts/hrom_census/quality_stratified_hypermode.py \
  --census "$CENSUS" \
  --genome-quality "$GENOME_QUALITY" \
  --species-metadata "$METADATA" \
  --outdir "$OUT"
