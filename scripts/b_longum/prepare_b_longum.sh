#!/bin/bash
# prepare_b_longum.sh — list B. longum genomes and create a manifest for
# downstream Syn2bANI all-vs-all comparison.
set -euo pipefail

GENOMES=${B_LONGUM_GENOMES:-/lustre1/g/aos_shihuang/Strain2b/data/JNU_genomes/genome2023/fna}
WORK=${B_LONGUM_WORK:-/lustre1/g/aos_shihuang/Syn2bANI-paper/results/b_longum_abfA}
mkdir -p "$WORK"

find "$GENOMES" -maxdepth 1 -name '*.fna' -o -name '*.fa' -o -name '*.fasta' | sort > "$WORK/genome_files.txt"
N=$(wc -l < "$WORK/genome_files.txt")
echo "Found $N genome files in $GENOMES"

# Create manifest: accession, file_path
awk '{n=$0; sub(/.*\//,"",n); sub(/\.(fna|fa|fasta)$/,"",n); print n"\t"$0}' "$WORK/genome_files.txt" > "$WORK/manifest.tsv"

# Create all-vs-all pair list (n*(n-1)/2)
awk 'NR==FNR{a[FNR]=$1;next}{for(i=1;i<FNR;i++)print a[i]"\t"$1}' "$WORK/genome_files.txt" "$WORK/genome_files.txt" > "$WORK/pairs_all_vs_all.tsv"
NPAIRS=$(wc -l < "$WORK/pairs_all_vs_all.tsv")
echo "Created $NPAIRS all-vs-all pairs"
