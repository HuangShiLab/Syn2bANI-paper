#!/bin/bash
#SBATCH --job-name=fetch_closed_inversions
#SBATCH --partition=amd
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/data/closed_inversions/logs/fetch_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/data/closed_inversions/logs/fetch_%j.err
# Task 7b from HPC_BRIEFING.md: fetch closed-genome cohort accessions.
set -euo pipefail
ROOT=/lustre1/g/aos_shihuang/Syn2bANI-paper
OUT=/lustre1/g/aos_shihuang/data/closed_inversions
mkdir -p "$OUT/logs"

bash "$ROOT/scripts/fetch_genomes_by_accession.sh" \
    "$OUT/accessions_to_download.txt" \
    "$OUT/genomes"
