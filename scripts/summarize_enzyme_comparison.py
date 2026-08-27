#!/usr/bin/env python3
import pandas as pd
import numpy as np

# Load mid-ANI matrix with FastANI reference
matrix = pd.read_csv("/lustre1/g/aos_shihuang/data/validation_mid_ani/mid_ani_matrix_gbrt_v4.tsv", sep="\t")
ref = matrix[['query', 'reference', 'fastani_ani', 'q_genus']].copy()

# Load enzyme comparison results
enz = pd.read_csv("/lustre1/g/aos_shihuang/data/validation_mid_ani/enzyme_comparison.tsv", sep="\t")

# Merge with FastANI reference
merged = enz.merge(ref, on=['query', 'reference'], how='inner')

# Filter out rows with errors
valid = merged[merged['raw_ani'].notna()].copy()
valid['error_raw'] = (valid['raw_ani'] - valid['fastani_ani']) * 100
valid['error_corrected'] = (valid['corrected_ani'] - valid['fastani_ani']) * 100

# Determine genus column after merge
genus_col = 'q_genus_y' if 'q_genus_y' in valid.columns else ('q_genus_x' if 'q_genus_x' in valid.columns else 'q_genus')

print("=" * 70)
print("Summary by enzyme mode (n = number of pairs)")
print("=" * 70)
summary = valid.groupby('enzyme_mode').agg(
    n=('query', 'count'),
    mean_raw_ani=('raw_ani', 'mean'),
    mean_shared_tags=('shared_tags', 'mean'),
    mean_af_q=('af_q', 'mean'),
    mean_af_r=('af_r', 'mean'),
    mean_fastani=('fastani_ani', 'mean'),
    mae_raw=('error_raw', lambda x: x.abs().mean()),
    mean_err_raw=('error_raw', 'mean'),
    mae_corrected=('error_corrected', lambda x: x.abs().mean()),
).reset_index()
print(summary.to_string(index=False, float_format='%.3f'))

print("\n" + "=" * 70)
print("By genus and enzyme mode")
print("=" * 70)
for genus in sorted(valid[genus_col].unique()):
    sub = valid[valid[genus_col] == genus]
    print(f"\n{genus}:")
    summ = sub.groupby('enzyme_mode').agg(
        n=('query', 'count'),
        mean_raw_ani=('raw_ani', 'mean'),
        mean_shared_tags=('shared_tags', 'mean'),
        mae_raw=('error_raw', lambda x: x.abs().mean()),
    ).reset_index()
    print(summ.to_string(index=False, float_format='%.3f'))
