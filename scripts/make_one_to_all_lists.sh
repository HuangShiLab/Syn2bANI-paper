#!/bin/bash
set -u
OUT=/Users/macstudio/Downloads/Syn2bANI-paper/results/efficiency_v8
REAL=/Users/macstudio/Downloads/Syn2bANI/prototype/realbench/genomes
DRAFT=/Users/macstudio/Downloads/Syn2bANI/prototype/draftbench/drafts
LISTS=$OUT/lists
mkdir -p "$LISTS"

POOL=(
  "$REAL/Ecoli_K12_MG1655.fasta"
  "$REAL/Ecoli_O157H7_Sakai.fasta"
  "$REAL/Ecoli_K12_W3110.fasta"
  "$REAL/Ecoli_CFT073.fasta"
  "$REAL/Ecoli_UTI89.fasta"
  "$REAL/Ecoli_BL21_DE3.fasta"
  "$REAL/Shigella_flexneri_301.fasta"
  "$REAL/Shigella_sonnei_Ss046.fasta"
  "$REAL/Escherichia_fergusonii.fasta"
  "$REAL/Salmonella_Typhimurium_LT2.fasta"
  "$REAL/Salmonella_Typhi_CT18.fasta"
  "$REAL/Klebsiella_pneumoniae_MGH78578.fasta"
  "$REAL/Enterobacter_cloacae_13047.fasta"
  "$REAL/Citrobacter_rodentium.fasta"
  "$DRAFT/GCA_001075925.fasta"
  "$DRAFT/GCA_001077875.fasta"
  "$DRAFT/GCA_001283205.fasta"
  "$DRAFT/GCA_001283245.fasta"
  "$DRAFT/GCA_001283605.fasta"
  "$DRAFT/GCA_001283865.fasta"
  "$DRAFT/GCA_001284145.fasta"
  "$DRAFT/GCA_001284645.fasta"
)

# one-to-all: first genome is query, next n are references
for n in 2 5 10 15 21; do
  FLIST=$LISTS/fasta_n${n}_one_to_all.txt
  # need query + n refs = n+1 genomes
  printf "%s\n" "${POOL[@]:0:$((n+1))}" > "$FLIST"
  echo "Wrote $FLIST with query + $n refs"
done
