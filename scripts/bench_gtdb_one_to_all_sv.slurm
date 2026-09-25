#!/bin/bash
#SBATCH --job-name=g1toall_sv
#SBATCH --partition=intel
#SBATCH --array=0-13
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/logs/g1toall_sv_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/logs/g1toall_sv_%A_%a.err

set -uo pipefail
BASE=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench
TOPHITS=$BASE/gtdb_one_to_all_tophits.tsv
QUERIES=$BASE/gtdb_one_to_all_queries.tsv
MODES=(syn2b_struct dnadiff)
N_QUERIES=$(($(wc -l < "$QUERIES") - 1))

IDX=$SLURM_ARRAY_TASK_ID
MODE=${MODES[$((IDX / N_QUERIES))]}
QID=$(tail -n +2 "$QUERIES" | sed -n "$((IDX % N_QUERIES + 1))p" | cut -f2)
LINE=$(grep -P "^${QID}\t" "$TOPHITS")
QFA=$(echo "$LINE" | cut -f2)
if [ "$MODE" = "syn2b_struct" ]; then
    REF=$(echo "$LINE" | cut -f3)
else
    REF=$(echo "$LINE" | cut -f4)
fi

ROWS=$BASE/results/rows_sv
mkdir -p "$ROWS"

echo "Task $IDX: mode=$MODE query=$QID ref=$REF"
python3 "$BASE/bench_gtdb_one_to_all.py" "$MODE" "$QFA" "$QID" "$ROWS" "$REF"
echo "Done task $IDX"
