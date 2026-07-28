#!/bin/bash -l
#SBATCH --job-name=s2b_v8_merge_p1
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/merge_v8_p1_%j.out
#SBATCH --error=logs/merge_v8_p1_%j.err

set -euo pipefail

source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani

OUT_DIR=/lustre1/g/aos_shihuang/data/gtdb-r207/matrix_v8_chunks
FINAL_MATRIX=results/matrix_gtdb_r207_100k_v8.tsv
FASTANI_PAIRS=results/pairs_fastani_subset_v8_100k.tsv
FASTANI_PAIRS_DIR=/lustre1/g/aos_shihuang/data/gtdb-r207/fastani_v8_pair_chunks
N_PER_BIN=200

export OUT_DIR FINAL_MATRIX FASTANI_PAIRS FASTANI_PAIRS_DIR N_PER_BIN

mkdir -p results "${FASTANI_PAIRS_DIR}"

python3 - <<'PY'
import pandas as pd
from pathlib import Path
import os
import sys

out_dir = Path(os.environ['OUT_DIR'])
final_matrix = os.environ['FINAL_MATRIX']
chunks = sorted(out_dir.glob('chunk_*.tsv'))
if not chunks:
    print('No chunk files found', file=sys.stderr)
    sys.exit(1)

dfs = []
for p in chunks:
    df = pd.read_csv(p, sep='\t', low_memory=False)
    dfs.append(df)
    print(f'Loaded {p}: {len(df)} rows')

combined = pd.concat(dfs, ignore_index=True)
combined = combined.drop_duplicates(subset=['query', 'reference'])
combined.to_csv(final_matrix, sep='\t', index=False, float_format='%.6f')
print(f'Merged matrix: {len(combined)} rows -> {final_matrix}')
PY

python3 scripts/sample_fastani_subset.py \
  --matrix "${FINAL_MATRIX}" \
  --output "${FASTANI_PAIRS}" \
  --n-per-bin "${N_PER_BIN}"

python3 scripts/split_pair_chunks.py \
  --pairs "${FASTANI_PAIRS}" \
  --output-dir "${FASTANI_PAIRS_DIR}" \
  --chunk-size 1000

echo "Phase 1 v8 merge complete. FastANI subset saved to ${FASTANI_PAIRS}"
