#!/bin/bash
#SBATCH --job-name=s2b_stra_full
#SBATCH --partition=amd
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=04:00:00
#SBATCH --output=strata_full_%j.out
#SBATCH --error=strata_full_%j.err

set -euo pipefail

cd /lustre1/g/aos_shihuang/Syn2bANI-paper

SYN2BANI=/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani
QL=results/gtdb_r207_100k_ql.txt
RL=results/gtdb_r207_100k_rl.txt
OUT=results/gtdb_r207_100k_strata.tsv

PANEL="BcgI,AlfI,AloI,FalI,BplI,Bsp24I,PpiI,PsrI,BsaXI,CjeI,CjePI"

echo "Starting full strata run at $(date) on $(hostname)"
echo "Pairs: $(wc -l < "$QL")"

time "$SYN2BANI" ani --ql "$QL" --rl "$RL" \
    -e "$PANEL" \
    --strata-out "$OUT" \
    -o results/gtdb_r207_100k_11enzyme.tsv \
    -p -t 64

echo "Finished at $(date)"
