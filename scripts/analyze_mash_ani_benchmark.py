#!/usr/bin/env python3
import pandas as pd
import numpy as np

VAL_DIR = '/lustre1/g/aos_shihuang/data/validation_mid_ani'

# Load FastANI reference
ref = pd.read_csv(f'{VAL_DIR}/mid_ani_matrix_gbrt_v4.tsv', sep='\t')[['query', 'reference', 'fastani_ani', 'skani_ani', 'q_genus']]

# Load new benchmark with mash_ani
mash = pd.read_csv(f'{VAL_DIR}/mid_ani_matrix_mash_ani.tsv', sep='\t')
mash = mash.merge(ref, on=['query', 'reference'], how='inner')
genus_col = 'q_genus_y' if 'q_genus_y' in mash.columns else ('q_genus_x' if 'q_genus_x' in mash.columns else 'q_genus')

for col in ['s2b_raw_ani', 's2b_mash_ani', 's2b_gbrt_ani', 'skani_ani']:
    err = (mash[col] - mash['fastani_ani']) * 100
    print(f"{col:20s}: MAE={err.abs().mean():.3f}%, mean_err={err.mean():.3f}%, rmse={np.sqrt((err**2).mean()):.3f}%")

print("\n=== Per-genus MAE ===")
for genus in sorted(mash[genus_col].unique()):
    sub = mash[mash[genus_col] == genus]
    print(f"\n{genus}:")
    for col in ['s2b_raw_ani', 's2b_mash_ani', 's2b_gbrt_ani', 'skani_ani']:
        err = (sub[col] - sub['fastani_ani']) * 100
        print(f"  {col:20s}: MAE={err.abs().mean():.3f}%")

print("\n=== Detailed pair comparison ===")
display = mash[['query', 'reference', genus_col, 'fastani_ani', 's2b_raw_ani', 's2b_mash_ani', 's2b_gbrt_ani', 'skani_ani']].copy()
for col in ['fastani_ani', 's2b_raw_ani', 's2b_mash_ani', 's2b_gbrt_ani', 'skani_ani']:
    display[col] = display[col] * 100
print(display.to_string(index=False, float_format='%.2f'))
