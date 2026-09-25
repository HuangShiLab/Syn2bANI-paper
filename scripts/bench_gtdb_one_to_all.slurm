#!/bin/bash
#SBATCH --job-name=g1toall_ani
#SBATCH --partition=intel
#SBATCH --array=0-34
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/logs/g1toall_ani_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/logs/g1toall_ani_%A_%a.err

set -uo pipefail
BASE=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench
QUERIES=$BASE/gtdb_one_to_all_queries.tsv
MODES=(syn2b_dist syn2b_search skani_dist skani_search fastani)
N_QUERIES=$(($(wc -l < "$QUERIES") - 1))

IDX=$SLURM_ARRAY_TASK_ID
MODE=${MODES[$((IDX / N_QUERIES))]}
LINE=$(tail -n +2 "$QUERIES" | sed -n "$((IDX % N_QUERIES + 1))p")
QFA=$(echo "$LINE" | cut -f6)
QID=$(echo "$LINE" | cut -f2)

ROWS=$BASE/results/rows_ani
mkdir -p "$ROWS"

echo "Task $IDX: mode=$MODE query=$QID fasta=$QFA"
python3 "$BASE/bench_gtdb_one_to_all.py" "$MODE" "$QFA" "$QID" "$ROWS"
echo "Done task $IDX"
