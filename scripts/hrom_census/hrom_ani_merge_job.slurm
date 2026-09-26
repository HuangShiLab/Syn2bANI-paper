#!/bin/bash
#SBATCH --job-name=hrom-merge-ani
#SBATCH --partition=amd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/databases/HROM/census/logs/merge_ani_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/databases/HROM/census/logs/merge_ani_%j.err

set -euo pipefail
REPO=/lustre1/g/aos_shihuang/Syn2b-paper
WORKDIR=/lustre1/g/aos_shihuang/databases/HROM/census
OUT=/lustre1/g/aos_shihuang/Syn2b-paper/results/census/hrom
mkdir -p "$WORKDIR/logs" "$OUT"
cd "$REPO"
python3 scripts/hrom_census/add_hrom_ani.py \
    --workdir "$WORKDIR" \
    --metadata "$WORKDIR/hrom_genome_metadata.tsv" \
    --structural "$OUT/hrom_within_species_sv.tsv.gz" \
    --out "$OUT/hrom_within_species_sv_ani.tsv.gz" \
    --summary "$OUT/ani_bin_summary.tsv" \
    --qc "$OUT/merge_ani_qc.md"
