#!/bin/bash -l
#SBATCH --job-name=s2b_mid_ani
#SBATCH --partition=amd
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/validate_mid_ani_%j.out
#SBATCH --error=logs/validate_mid_ani_%j.err

set -euo pipefail

source /group/aos_shihuang/conda/etc/profile.d/conda.sh
conda activate syn2bani

VAL_DIR=/lustre1/g/aos_shihuang/data/validation_mid_ani
GENOMES_DIR="${VAL_DIR}/genomes"
PAIRS_FILE="${VAL_DIR}/mid_ani_pairs.tsv"
FASTANI_MATRIX="${VAL_DIR}/mid_ani_matrix_fastani.tsv"
FILTERED_PAIRS="${VAL_DIR}/mid_ani_pairs_85_95.tsv"
TOOLS_MATRIX="${VAL_DIR}/mid_ani_matrix_tools.tsv"
FINAL_MATRIX="${VAL_DIR}/mid_ani_matrix.tsv"
SYN2BANI=/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani
MODEL=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gbrt_model_v4_100k.pkl
REPORT=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/VALIDATION_MID_ANI_REPORT.md
PLOT_DIR=/lustre1/g/aos_shihuang/Syn2bANI-paper/figures/validation_mid_ani

mkdir -p "${VAL_DIR}" logs "${PLOT_DIR}"
cd /lustre1/g/aos_shihuang/Syn2bANI-paper

echo "=== Step 1: Build within-genus pairs ==="
python3 scripts/build_mid_ani_pairs.py \
  --manifest "${VAL_DIR}/manifest.tsv" \
  --output "${PAIRS_FILE}"

echo ""
echo "=== Step 2: Run FastANI on all within-genus pairs ==="
/usr/bin/time -v python3 scripts/run_benchmark_matrix_v2.py \
  --pairs "${PAIRS_FILE}" \
  --genomes "${GENOMES_DIR}" \
  --syn2bani "${SYN2BANI}" \
  --skani skani \
  --fastani fastANI \
  --output "${FASTANI_MATRIX}" \
  --threads 8 \
  --tools fastani \
  --chunk-size 100

echo ""
echo "=== Step 3: Filter pairs with 0.85 <= FastANI ANI <= 0.95 ==="
python3 - <<'PY'
import pandas as pd
from pathlib import Path

val_dir = Path('/lustre1/g/aos_shihuang/data/validation_mid_ani')
df = pd.read_csv(val_dir / 'mid_ani_matrix_fastani.tsv', sep='\t')
print(f'Pairs before filtering: {len(df)}')
print(f'FastANI ANI range: {df["fastani_ani"].min():.4f} - {df["fastani_ani"].max():.4f}')

mask = df['fastani_ani'].notna() & (df['fastani_ani'] >= 0.85) & (df['fastani_ani'] <= 0.95)
filtered = df[mask].copy()
print(f'Pairs in 85-95% ANI range: {len(filtered)}')

# Keep only pair-identification columns for the next tool run
keep_cols = ['query', 'reference', 'label', 'q_species', 'r_species',
             'q_genus', 'r_genus', 'q_category', 'r_category']
filtered = filtered[[c for c in keep_cols if c in filtered.columns]]
filtered.to_csv(val_dir / 'mid_ani_pairs_85_95.tsv', sep='\t', index=False)
PY

echo ""
echo "=== Step 4: Run skani + Syn2bANI on filtered pairs ==="
/usr/bin/time -v python3 scripts/run_benchmark_matrix_v2.py \
  --pairs "${FILTERED_PAIRS}" \
  --genomes "${GENOMES_DIR}" \
  --syn2bani "${SYN2BANI}" \
  --skani skani \
  --fastani fastANI \
  --output "${TOOLS_MATRIX}" \
  --threads 8 \
  --tools skani,syn2bani \
  --chunk-size 100

echo ""
echo "=== Step 5: Merge FastANI reference with tool predictions ==="
python3 - <<'PY'
import pandas as pd
from pathlib import Path

val_dir = Path('/lustre1/g/aos_shihuang/data/validation_mid_ani')
fastani = pd.read_csv(val_dir / 'mid_ani_matrix_fastani.tsv', sep='\t')
tools = pd.read_csv(val_dir / 'mid_ani_matrix_tools.tsv', sep='\t')

# Keep FastANI result plus pair metadata from the full matrix; merge tool preds
keep = ['query', 'reference', 'label', 'q_species', 'r_species',
        'q_genus', 'r_genus', 'q_category', 'r_category', 'fastani_ani']
fastani = fastani[[c for c in keep if c in fastani.columns]]

tool_cols = ['query', 'reference', 'skani_ani', 'skani_align_frac',
             's2b_raw_ani', 's2b_gbrt_ani', 's2b_shared_tags',
             's2b_af_q', 's2b_af_r', 's2b_ref_gc']
tools = tools[[c for c in tool_cols if c in tools.columns]]

merged = fastani.merge(tools, on=['query', 'reference'], how='inner')
merged.to_csv(val_dir / 'mid_ani_matrix.tsv', sep='\t', index=False, float_format='%.6f')
print(f'Merged matrix saved: {len(merged)} rows')
PY

echo ""
echo "=== Step 6: Evaluate GBRT v4 ==="
python3 scripts/evaluate_validation.py \
  --matrix "${FINAL_MATRIX}" \
  --model "${MODEL}" \
  --output "${FINAL_MATRIX%.tsv}_gbrt_v4.tsv" \
  --report "${REPORT}"

echo ""
echo "=== Step 7: Plot validation results ==="
python3 scripts/plot_validation.py \
  --matrix "${FINAL_MATRIX%.tsv}_gbrt_v4.tsv" \
  --out-dir "${PLOT_DIR}"

echo ""
echo "=== Validation complete ==="
ls -lh "${FINAL_MATRIX}"
ls -lh "${REPORT}"
ls -lh "${PLOT_DIR}"
