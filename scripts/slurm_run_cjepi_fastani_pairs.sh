#!/bin/bash -l
#SBATCH --job-name=s2b_cjepi
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=logs/run_cjepi_fastani_pairs_%j.out
#SBATCH --error=logs/run_cjepi_fastani_pairs_%j.err

set -euo pipefail

source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani

cd /lustre1/g/aos_shihuang/Syn2bANI-paper

python3 scripts/run_benchmark_enzyme.py \
  --pairs /lustre1/g/aos_shihuang/data/gtdb-r207/fastani_pairs_728.tsv \
  --genomes /lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all \
  --syn2bani /lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani \
  --enzyme CjePI \
  --output /lustre1/g/aos_shihuang/data/gtdb-r207/fastani_pairs_728_cjepi.tsv \
  --threads 16

echo "=== CjePI benchmark complete ==="
ls -lh /lustre1/g/aos_shihuang/data/gtdb-r207/fastani_pairs_728_cjepi.tsv
