#!/bin/bash -l
#SBATCH --job-name=s2b_dbscale
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --time=8:00:00
#SBATCH --output=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/db_scale/logs/slurm_%j.out
#SBATCH --error=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/db_scale/logs/slurm_%j.err

set -euo pipefail
export THREADS=32
DIR=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/db_scale/scripts
source "$DIR/common.sh"

echo "host=$(hostname) job=$SLURM_JOB_ID start=$(date -u +%FT%TZ)"
lscpu | grep -E 'Model name|^CPU\(s\)' || true

bash "$DIR/bench_sketch.sh" 2000 2
bash "$DIR/bench_sketch.sh" 5000 1
bash "$DIR/bench_triangle.sh" 2000 1
bash "$DIR/bench_search.sh" 2
bash "$DIR/bench_accuracy.sh" 500
echo "ALL DONE $(date -u +%FT%TZ)"
