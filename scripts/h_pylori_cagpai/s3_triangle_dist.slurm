#!/bin/bash
# s3_triangle_dist.slurm — all-vs-all ANI + af for the cohort (edge list).
#SBATCH --job-name=hp_triangle
#SBATCH --output=hp_triangle_%j.out
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00

# Env: GENOMES (dir with *.fna), OUT (output tsv), S2B (syn2bani >= 4d679cf)
GENOMES=${GENOMES:-$PWD/genomes}
OUT=${OUT:-$PWD/hp_triangle.tsv}
S2B=${S2B:-$TOOLS/syn2bani/syn2bani}

ls "$GENOMES"/*.fna > genomes_list.txt
# --verbose adds af_query/af_reference/synteny columns; --edge-list gives
# ani-style rows for every refined pair.
"$S2B" triangle --ql genomes_list.txt --edge-list --verbose -p -o "$OUT"
echo "pairs: $(($(wc -l < "$OUT") - 1))"
