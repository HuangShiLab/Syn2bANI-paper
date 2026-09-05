#!/bin/bash
#SBATCH --job-name=s2b_struct_bench
#SBATCH --partition=intel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench/logs/%A_%a.log

set -u
BASE=/lustre1/g/aos_shihuang/Syn2bANI-paper-bench
OUT=$BASE/results
LOGS=$BASE/logs
mkdir -p $OUT $LOGS

PY=$BASE/bench_syn2b_struct.py

NS=(2 5 10 15 22)
IDX=$((SLURM_ARRAY_TASK_ID - 1))
N_IDX=$((IDX / 3))
REP_IDX=$((IDX % 3))
N=${NS[$N_IDX]}
REP=$((REP_IDX + 1))

FLIST=$BASE/lists/fasta_n$N.txt
if [ ! -f "$FLIST" ]; then
    echo "List not found: $FLIST" >&2
    exit 1
fi

RES=$OUT/syn2b_struct_benchmark.tsv
if [ ! -f "$RES" ]; then
    printf "mode\tn_genomes\tn_pairs\trep\tstruct_wall_s\tn_ok\n" > "$RES"
fi

/usr/bin/time -v python3 "$PY" "$FLIST" "$N" "$REP" "$RES" \
    > "$LOGS/syn2b_struct_n${N}_r${REP}.out" 2> "$LOGS/syn2b_struct_n${N}_r${REP}.err"

echo "Done task $SLURM_ARRAY_TASK_ID: n=$N rep=$REP"
