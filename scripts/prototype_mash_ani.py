#!/usr/bin/env python3
"""Prototype Mash-like ANI estimators using existing enzyme comparison data."""
import pandas as pd
import numpy as np
from pathlib import Path

VAL_DIR = Path('/lustre1/g/aos_shihuang/data/validation_mid_ani')

# Load reference FastANI
ref = pd.read_csv(VAL_DIR / 'mid_ani_matrix_gbrt_v4.tsv', sep='\t')[['query', 'reference', 'fastani_ani', 'q_genus']]

# Load single-enzyme comparison
single = pd.read_csv(VAL_DIR / 'enzyme_comparison.tsv', sep='\t')
single = single[single['raw_ani'].notna()].copy()
single = single.merge(ref, on=['query', 'reference'], how='inner')

# Load multi-enzyme results
multi = pd.read_csv(VAL_DIR / 'multi_enzyme_mid_ani.tsv', sep='\t')
multi = multi[multi['raw_ani'].notna()].copy()
multi = multi.merge(ref, on=['query', 'reference'], how='inner')
multi['enzyme_mode'] = 'multi'

# Combine
combined = pd.concat([single, multi], ignore_index=True)

# Effective tag lengths for each enzyme (from registry.rs)
TAG_LENGTHS = {
    'BcgI': 32,
    'BsaXI': 27,
    'CjeI': 28,
    'CjePI': 27,
    'FalI': 27,
    'HaeIV': 27,
    'multi': 28,  # approximate average across 16 enzymes
}


def compute_mash_ani(af_q, af_r, tag_len):
    """Mash-like ANI from bidirectional containment."""
    af_geo = np.sqrt(np.maximum(af_q, 1e-10) * np.maximum(af_r, 1e-10))
    return 1.0 + np.log(af_geo) / tag_len


def compute_mash_ani_jaccard(af_q, af_r, tag_len):
    """Estimate Jaccard from containment and convert to Mash distance."""
    # Jaccard ≈ C_q * C_r / (C_q + C_r - C_q*C_r)
    c1 = np.maximum(af_q, 1e-10)
    c2 = np.maximum(af_r, 1e-10)
    jaccard = (c1 * c2) / (c1 + c2 - c1 * c2)
    # Mash distance: D = -ln(2J/(1+J)) / k
    mash_dist = -np.log(2 * jaccard / (1 + jaccard)) / tag_len
    return np.maximum(0.0, 1.0 - mash_dist)


# Apply formulas
for name, tag_len in TAG_LENGTHS.items():
    mask = combined['enzyme_mode'] == name
    if mask.sum() == 0:
        continue
    combined.loc[mask, 'mash_ani'] = compute_mash_ani(
        combined.loc[mask, 'af_q'], combined.loc[mask, 'af_r'], tag_len
    )
    combined.loc[mask, 'mash_ani_jaccard'] = compute_mash_ani_jaccard(
        combined.loc[mask, 'af_q'], combined.loc[mask, 'af_r'], tag_len
    )

# Compute errors
combined['error_raw'] = (combined['raw_ani'] - combined['fastani_ani']) * 100
combined['error_corrected'] = (combined['corrected_ani'] - combined['fastani_ani']) * 100
combined['error_mash'] = (combined['mash_ani'] - combined['fastani_ani']) * 100
combined['error_mash_jaccard'] = (combined['mash_ani_jaccard'] - combined['fastani_ani']) * 100

print("=" * 90)
print("MAE by enzyme mode and estimator")
print("=" * 90)
summary = combined.groupby('enzyme_mode').agg(
    n=('query', 'count'),
    mean_fastani=('fastani_ani', 'mean'),
    mean_raw=('raw_ani', 'mean'),
    mae_raw=('error_raw', lambda x: x.abs().mean()),
    mae_corrected=('error_corrected', lambda x: x.abs().mean()),
    mae_mash=('error_mash', lambda x: x.abs().mean()),
    mae_mash_jaccard=('error_mash_jaccard', lambda x: x.abs().mean()),
).reset_index()
print(summary.to_string(index=False, float_format='%.3f'))

print("\n" + "=" * 90)
print("Detailed: multi-enzyme results")
print("=" * 90)
multi_sub = combined[combined['enzyme_mode'] == 'multi'][['query', 'reference', 'fastani_ani', 'raw_ani', 'mash_ani', 'mash_ani_jaccard', 'corrected_ani', 'af_q', 'af_r', 'shared_tags']]
print(multi_sub.to_string(index=False, float_format='%.4f'))

print("\n" + "=" * 90)
print("Best single enzyme (lowest Mash ANI MAE)")
print("=" * 90)
single_sub = combined[combined['enzyme_mode'] != 'multi']
best = single_sub.groupby('enzyme_mode')['error_mash'].apply(lambda x: x.abs().mean()).sort_values()
print(best.to_string(float_format='%.3f'))
