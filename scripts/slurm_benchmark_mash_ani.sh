#!/bin/bash -l
#SBATCH --job-name=s2b_mash_ani
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=logs/benchmark_mash_ani_%j.out
#SBATCH --error=logs/benchmark_mash_ani_%j.err

set -euo pipefail

source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani

cd /lustre1/g/aos_shihuang/Syn2bANI-paper

VAL_DIR=/lustre1/g/aos_shihuang/data/validation_mid_ani
SYN2BANI=/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani

echo "=== Benchmark Mash-like ANI on mid-ANI pairs ==="
/usr/bin/time -v python3 scripts/run_benchmark_matrix_v2.py \
  --pairs "${VAL_DIR}/mid_ani_pairs_85_95.tsv" \
  --genomes "${VAL_DIR}/genomes" \
  --syn2bani "${SYN2BANI}" \
  --skani skani \
  --fastani fastANI \
  --output "${VAL_DIR}/mid_ani_matrix_mash_ani.tsv" \
  --threads 8 \
  --tools syn2bani \
  --chunk-size 100

echo "=== Benchmark complete ==="
ls -lh "${VAL_DIR}/mid_ani_matrix_mash_ani.tsv"
