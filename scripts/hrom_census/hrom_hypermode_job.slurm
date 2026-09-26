#!/bin/bash
#SBATCH --job-name=hrom-hypermodes
#SBATCH --partition=amd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/databases/HROM/census/logs/hypermodes_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/databases/HROM/census/logs/hypermodes_%j.err

set -euo pipefail
REPO=/lustre1/g/aos_shihuang/Syn2b-paper
CENSUS=/lustre1/g/aos_shihuang/Syn2b-paper/results/census/hrom
METADATA=/lustre1/g/aos_shihuang/databases/HROM/data/genome_catalog/HROM_representative_genome_metadata.tsv
OUT=$CENSUS/hypermodes
mkdir -p "$OUT"
cd "$REPO"
python3 scripts/hrom_census/analyze_hrom_hypermodes.py \
    --census "$CENSUS/hrom_within_species_sv_ani.tsv.gz" \
    --species-metadata "$METADATA" \
    --outdir "$OUT"
python3 scripts/hrom_census/plot_hrom_hypermodes.py \
    --points "$OUT/hypermode_figure_points.tsv.gz" \
    --enrichment "$OUT/hypermode_species_enrichment.tsv" \
    --report "$OUT/hypermode_analysis_report.json" \
    --outstem "$OUT/fig_hrom_hypermodes"
