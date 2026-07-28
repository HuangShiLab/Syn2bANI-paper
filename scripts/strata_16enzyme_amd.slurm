#!/bin/bash
#SBATCH --job-name=s2b_strata_amd
#SBATCH --partition=amd
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=24:00:00
#SBATCH --output=strata_amd_%j.out
#SBATCH --error=strata_amd_%j.err

set -euo pipefail

cd /lustre1/g/aos_shihuang/Syn2bANI-paper

SYN2BANI=/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani
QL=results/gtdb_r207_100k_ql.txt
RL=results/gtdb_r207_100k_rl.txt
OUT=results/gtdb_r207_100k_strata.tsv

PANEL="BcgI,AlfI,AloI,FalI,BplI,Bsp24I,PpiI,PsrI,BsaXI,CjeI,CjePI"

echo "Starting strata run at $(date)"
echo "Pairs: $(wc -l < "$QL")"

$SYN2BANI ani --ql "$QL" --rl "$RL" \
    -e "$PANEL" \
    --strata-out "$OUT" \
    -o results/gtdb_r207_100k_11enzyme.tsv \
    -p -t 32

echo "Finished at $(date)"
