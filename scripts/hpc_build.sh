#!/bin/bash -l
#SBATCH --job-name=s2b_build
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=logs/build_%j.out
#SBATCH --error=logs/build_%j.err

set -euo pipefail

JOB_ID="${SLURM_JOB_ID:-login}"
LOG_SUFFIX="${JOB_ID}"

echo "Building Syn2bANI on $(hostname) at $(date)"

cd /lustre1/g/aos_shihuang/Syn2bANI
mkdir -p logs

cargo build --release 2>&1 | tee -a logs/build_${LOG_SUFFIX}.log
cargo test --release 2>&1 | tee -a logs/build_${LOG_SUFFIX}.log

echo "Build complete at $(date)"
ls -lh target/release/syn2bani | tee -a logs/build_${LOG_SUFFIX}.log
