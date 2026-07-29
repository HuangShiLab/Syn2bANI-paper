#!/bin/bash
#SBATCH --job-name=s2b_stra_exc
#SBATCH --partition=amd
#SBATCH --array=1-10
#SBATCH --exclusive
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=04:00:00
#SBATCH --output=strata_exc_%A_%a.out
#SBATCH --error=strata_exc_%A_%a.err

set -euo pipefail

cd /lustre1/g/aos_shihuang/Syn2bANI-paper

SYN2BANI=/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani
IDX=$(printf "%02d" "$SLURM_ARRAY_TASK_ID")
QL="results/strata_chunks/ql_${IDX}.txt"
RL="results/strata_chunks/rl_${IDX}.txt"
OUT="results/strata_chunks/strata_${IDX}.tsv"
ANI_OUT="results/strata_chunks/ani_${IDX}.tsv"

PANEL="BcgI,AlfI,AloI,FalI,BplI,Bsp24I,PpiI,PsrI,BsaXI,CjeI,CjePI"

echo "Starting chunk $IDX at $(date) on $(hostname)"
echo "Pairs: $(wc -l < "$QL")"

time "$SYN2BANI" ani --ql "$QL" --rl "$RL" \
    -e "$PANEL" \
    --strata-out "$OUT" \
    -o "$ANI_OUT" \
    -p -t 64

echo "Finished chunk $IDX at $(date)"
