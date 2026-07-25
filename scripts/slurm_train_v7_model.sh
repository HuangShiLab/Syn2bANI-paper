#!/bin/bash -l
#SBATCH --job-name=s2b_v7_train
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=logs/s2b_v7_train_%j.out
#SBATCH --error=logs/s2b_v7_train_%j.err
#SBATCH --dependency=afterany:?FEATURES_JOBID?

set -euo pipefail

source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani

GTDB_DIR=/lustre1/g/aos_shihuang/databases/GTDB/GTDBr207
FEATURES_DIR="${GTDB_DIR}/v7_features"
MATRIX="${GTDB_DIR}/matrix_gtdb_r207_100k.tsv"
PAPER_DIR=/lustre1/g/aos_shihuang/Syn2bANI-paper
OUTPUT_DIR="${PAPER_DIR}/results/gbrt_v7"

mkdir -p "${OUTPUT_DIR}" logs

echo "=== Merging v7 feature chunks with FastANI reference ==="
python3 - <<PY
import sys
import pandas as pd
from pathlib import Path

features_dir = Path('${FEATURES_DIR}')
matrix_path = Path('${MATRIX}')
output_path = Path('${OUTPUT_DIR}/gtdb_r207_v7_training_matrix.tsv')

chunks = sorted(features_dir.glob('features_chunk_*.tsv'))
print(f'Found {len(chunks)} feature chunks', file=sys.stderr)

feat_list = []
for chunk in chunks:
    df = pd.read_csv(chunk, sep='\t', low_memory=False)
    feat_list.append(df)
features = pd.concat(feat_list, ignore_index=True)
print(f'Features rows: {len(features)}', file=sys.stderr)

matrix = pd.read_csv(matrix_path, sep='\t', low_memory=False)
# Keep only necessary columns from matrix
keep_cols = ['query', 'reference', 'label', 'q_genus', 'r_genus', 'q_species', 'r_species',
             'skani_ani', 'skani_align_frac', 'fastani_ani']
keep_cols = [c for c in keep_cols if c in matrix.columns]
matrix = matrix[keep_cols]

merged = matrix.merge(features, on=['query', 'reference'], how='inner')
print(f'Merged rows: {len(merged)}', file=sys.stderr)

# Filter to rows with valid FastANI reference
merged = merged[merged['fastani_ani'].notna() & (merged['fastani_ani'] > 0)]
print(f'Rows with valid FastANI: {len(merged)}', file=sys.stderr)

output_path.parent.mkdir(parents=True, exist_ok=True)
merged.to_csv(output_path, sep='\t', index=False, float_format='%.6f')
print(f'Wrote training matrix to {output_path}', file=sys.stderr)
PY

echo ""
echo "=== Training GBRT v7 model ==="
python3 "${PAPER_DIR}/scripts/train_gbrt_v7.py" \
  --matrix "${OUTPUT_DIR}/gtdb_r207_v7_training_matrix.tsv" \
  --output "${OUTPUT_DIR}/gbrt_model_v7.json" \
  --report "${OUTPUT_DIR}/gbrt_v7_report.txt"

echo ""
echo "=== Training complete ==="
ls -lh "${OUTPUT_DIR}"
