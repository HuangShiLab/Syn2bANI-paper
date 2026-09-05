#!/bin/bash
#SBATCH --job-name=gtdb_1toall
#SBATCH --partition=intel
#SBATCH --array=0-34
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/logs/gtdb_1toall_%A_%a.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/logs/gtdb_1toall_%A_%a.err

set -euo pipefail
BASE=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench
OUT=$BASE/results
LOGS=$BASE/logs
mkdir -p $OUT $LOGS

QUERIES=$BASE/gtdb_one_to_all_queries.tsv
MODES=(syn2b_fasta syn2b_search skani_fasta skani_sketch fastani)
N_MODES=${#MODES[@]}
N_QUERIES=$(tail -n +2 "$QUERIES" | wc -l)

IDX=$SLURM_ARRAY_TASK_ID
MODE_IDX=$((IDX / N_QUERIES))
QUERY_IDX=$((IDX % N_QUERIES))
MODE=${MODES[$MODE_IDX]}
QUERY_LINE=$(tail -n +2 "$QUERIES" | sed -n "$((QUERY_IDX + 1))p")
QUERY_FA=$(echo "$QUERY_LINE" | cut -f6)
QUERY_NAME=$(echo "$QUERY_LINE" | cut -f2)

RES=$OUT/gtdb_one_to_all_scaling.tsv
if [ ! -f "$RES" ]; then
    printf "mode\tquery\twall_s\tpeak_rss_mb\n" > "$RES"
fi

echo "Task $IDX: mode=$MODE query=$QUERY_NAME"
/usr/bin/time -v python3 "$BASE/bench_gtdb_one_to_all.py" "$MODE" "$QUERY_FA" "$RES" \
    > "$LOGS/gtdb_1toall_${MODE}_${QUERY_NAME}.out" 2> "$LOGS/gtdb_1toall_${MODE}_${QUERY_NAME}.err"

echo "Done task $IDX"
