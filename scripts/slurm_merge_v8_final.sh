#!/bin/bash -l
#SBATCH --job-name=s2b_v8_merge_final
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/merge_v8_final_%j.out
#SBATCH --error=logs/merge_v8_final_%j.err

set -euo pipefail

source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani

MATRIX=results/matrix_gtdb_r207_100k_v8.tsv
FASTANI_DIR=/lustre1/g/aos_shihuang/data/gtdb-r207/fastani_v8_chunks
FINAL_MATRIX=results/matrix_gtdb_r207_100k_v8_final.tsv

export MATRIX FASTANI_DIR FINAL_MATRIX

python3 - <<'PY'
import pandas as pd
from pathlib import Path
import os
import sys

matrix_path = os.environ['MATRIX']
fastani_dir = Path(os.environ['FASTANI_DIR'])
final_path = os.environ['FINAL_MATRIX']

matrix = pd.read_csv(matrix_path, sep='\t', low_memory=False)
print(f'Loaded matrix: {len(matrix)} rows')

chunks = sorted(fastani_dir.glob('chunk_*.tsv'))
if not chunks:
    print('No FastANI chunk files found', file=sys.stderr)
    matrix.to_csv(final_path, sep='\t', index=False, float_format='%.6f')
    print(f'No FastANI results to merge. Saved {final_path}')
    sys.exit(0)

fastani_dfs = []
for p in chunks:
    df = pd.read_csv(p, sep='\t', low_memory=False)
    fastani_dfs.append(df)
    print(f'Loaded FastANI chunk {p}: {len(df)} rows')

fastani = pd.concat(fastani_dfs, ignore_index=True)
fastani = fastani.drop_duplicates(subset=['query', 'reference'])
fastani = fastani[['query', 'reference', 'fastani_ani']]

merged = matrix.merge(fastani, on=['query', 'reference'], how='left', suffixes=('', '_fa'))
# If matrix already has fastani_ani from a previous run, prefer the new one
if 'fastani_ani_fa' in merged.columns:
    merged['fastani_ani'] = merged['fastani_ani_fa'].combine_first(merged.get('fastani_ani'))
    merged = merged.drop(columns=['fastani_ani_fa'])

merged.to_csv(final_path, sep='\t', index=False, float_format='%.6f')
print(f'Merged final matrix: {len(merged)} rows -> {final_path}')
print(f'Rows with FastANI reference: {merged["fastani_ani"].notna().sum()}')
PY

echo "Final v8 merge complete: ${FINAL_MATRIX}"
