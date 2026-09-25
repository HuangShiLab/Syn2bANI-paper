#!/bin/bash
#SBATCH --job-name=g1toall_loop
#SBATCH --partition=intel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/logs/g1toall_loop_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/logs/g1toall_loop_%j.err

# Submit with: sbatch --export=ALL,MODE=syn2b_dist,QSTART=0,QEND=7 bench_gtdb_one_to_all_loop.slurm
# QSTART inclusive, QEND exclusive (row indices into the queries file, header excluded).

set -uo pipefail
BASE=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench
QUERIES=$BASE/gtdb_one_to_all_queries.tsv
MODE=${MODE:?MODE required}
QSTART=${QSTART:-0}
QEND=${QEND:-7}

ROWS=$BASE/results/rows_ani
mkdir -p "$ROWS"

mapfile -t LINES < <(tail -n +2 "$QUERIES")
for ((i=QSTART; i<QEND; i++)); do
    LINE="${LINES[$i]}"
    QFA=$(echo "$LINE" | cut -f6)
    QID=$(echo "$LINE" | cut -f2)
    echo "=== mode=$MODE query=$QID ($i) ==="
    python3 "$BASE/bench_gtdb_one_to_all.py" "$MODE" "$QFA" "$QID" "$ROWS" || echo "FAILED $MODE $QID rc=$?"
done
echo "Loop done: $MODE [$QSTART,$QEND)"
