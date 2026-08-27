#!/usr/bin/env python3
"""Evaluate a GBRT v5 pickle model on a validation matrix."""
import argparse
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, mean_squared_error


def compute_mash_ani(af_q, af_r, tag_len=32.0):
    c1 = np.maximum(af_q, 1e-10)
    c2 = np.maximum(af_r, 1e-10)
    containment_geo = np.sqrt(c1 * c2)
    return np.clip(1.0 + np.log(containment_geo) / tag_len, 0.0, 1.0)


def build_features(df):
    df = df[df['s2b_raw_ani'].notna()].copy()
    if 's2b_mash_ani' not in df.columns:
        df['s2b_mash_ani'] = compute_mash_ani(df['s2b_af_q'].values, df['s2b_af_r'].values)

    df['raw_ani'] = df['s2b_raw_ani']
    df['mash_ani'] = df['s2b_mash_ani']
    df['shared_log'] = np.log1p(df['s2b_shared_tags'].fillna(0).clip(lower=0))
    df['af_q'] = df['s2b_af_q'].fillna(df['s2b_af_q'].median()).clip(0, 1)
    df['af_r'] = df['s2b_af_r'].fillna(df['s2b_af_r'].median()).clip(0, 1)

    feature_cols = ['raw_ani', 'mash_ani', 'shared_log', 'af_q', 'af_r']
    X_df = df[feature_cols].fillna(df[feature_cols].median())
    return X_df.values, df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', required=True)
    parser.add_argument('--matrix', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--label', default='validation')
    args = parser.parse_args()

    with open(args.model, 'rb') as f:
        model = pickle.load(f)

    df = pd.read_csv(args.matrix, sep='\t', low_memory=False)
    X, df = build_features(df)

    bias = model.predict(X)
    df['gbrt_v5_bias'] = bias
    df['gbrt_v5_ani'] = (df['s2b_raw_ani'] + bias).clip(0, 1)

    ref_col = 'fastani_ani'
    if ref_col in df.columns:
        mask = df[ref_col].notna() & df['gbrt_v5_ani'].notna()
        sub = df[mask]

        methods = {
            'Syn2bANI raw': 's2b_raw_ani',
            'Syn2bANI mash': 's2b_mash_ani',
            'Syn2bANI GBRT v5': 'gbrt_v5_ani',
        }
        if 's2b_gbrt_ani' in df.columns:
            methods['Syn2bANI GBRT v4'] = 's2b_gbrt_ani'
        if 'skani_ani' in df.columns:
            methods['skani'] = 'skani_ani'

        print(f"\n=== {args.label} Metrics (n={len(sub)}) ===")
        print(f"{'Method':<25} {'MAE':>10} {'RMSE':>10}")
        for name, col in methods.items():
            if col in sub.columns:
                m = sub[col]
                r = sub[ref_col]
                mae = mean_absolute_error(r, m) * 100
                rmse = np.sqrt(mean_squared_error(r, m)) * 100
                print(f"{name:<25} {mae:>10.3f}% {rmse:>10.3f}%")

    df.to_csv(args.output, sep='\t', index=False, float_format='%.6f')
    print(f"\nSaved: {args.output}")


if __name__ == '__main__':
    main()
