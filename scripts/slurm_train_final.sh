#!/bin/bash -l
#SBATCH --job-name=s2b_train
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/train_%j.out
#SBATCH --error=logs/train_%j.err

set -euo pipefail

# Activate conda environment
source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani

# Configuration
MATRIX_DIR=/lustre1/g/aos_shihuang/data/gtdb-r207/matrix_chunks
FASTANI_DIR=/lustre1/g/aos_shihuang/data/gtdb-r207/fastani_chunks
SKANI_MATRIX=results/matrix_gtdb_r207_100k_skani.tsv
FINAL_MATRIX=results/matrix_gtdb_r207_100k.tsv
MODEL_JSON=results/gbrt_model_v7_100k.json
MODEL_REPORT=results/gbrt_v7_100k_report.txt
EVAL_JSON=results/evaluation_gtdb_r207_100k.json
FIG_DIR=figures

export FASTANI_DIR SKANI_MATRIX FINAL_MATRIX

mkdir -p results "${FIG_DIR}"

# Merge FastANI chunk results
python3 - <<'PY'
import pandas as pd
from pathlib import Path
import os
import sys

fastani_dir = Path(os.environ['FASTANI_DIR'])
skani_matrix = os.environ['SKANI_MATRIX']
final_matrix = os.environ['FINAL_MATRIX']

chunks = sorted(fastani_dir.glob('chunk_*.tsv'))
if not chunks:
    print('No FastANI chunk files found', file=sys.stderr)
    sys.exit(1)

fastani_parts = []
for p in chunks:
    df = pd.read_csv(p, sep='\t', low_memory=False)
    fastani_parts.append(df[['query', 'reference', 'fastani_ani']])
    print(f'Loaded FastANI chunk {p}: {len(df)} rows')

fastani_df = pd.concat(fastani_parts, ignore_index=True)
fastani_df = fastani_df.drop_duplicates(subset=['query', 'reference'])
print(f'Total FastANI rows: {len(fastani_df)}')

# Load the Syn2bANI + skani matrix and merge FastANI reference
skani_df = pd.read_csv(skani_matrix, sep='\t', low_memory=False)
if 'fastani_ani' in skani_df.columns:
    skani_df = skani_df.drop(columns=['fastani_ani'])
merged = skani_df.merge(fastani_df, on=['query', 'reference'], how='left')
merged.to_csv(final_matrix, sep='\t', index=False, float_format='%.6f')
print(f'Final matrix: {len(merged)} rows, FastANI valid: {merged["fastani_ani"].notna().sum()}')
PY

# Train clean GBRT v7 model on the FastANI reference
python3 scripts/train_gbrt_v7.py \
  --matrix "${FINAL_MATRIX}" \
  --output "${MODEL_JSON}" \
  --report "${MODEL_REPORT}"

# Evaluate
python3 scripts/evaluate_gtdb_r207.py \
  --matrix "${FINAL_MATRIX}" \
  --output "${EVAL_JSON}" \
  --figures "${FIG_DIR}" \
  --model "${MODEL_JSON%.json}.pkl"

echo "Final merge, training, and evaluation complete"
