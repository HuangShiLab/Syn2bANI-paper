#!/usr/bin/env python3
"""Evaluate GBRT v5 models for BcgI, CjePI, or combined features on a validation matrix."""
import argparse
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, mean_squared_error


def build_features(df, mode):
    """Build feature matrix X and metadata."""
    prefixes = []
    if mode in ('bcgi', 'combined'):
        prefixes.append('bcgi')
    if mode in ('cjepi', 'combined'):
        prefixes.append('cjepi')

    base_prefix = prefixes[0]
    df['base_raw_ani'] = df[f'{base_prefix}_raw_ani']

    feature_cols = []
    for prefix in prefixes:
        df[f'{prefix}_shared_log'] = np.log1p(df[f'{prefix}_shared_tags'].fillna(0).clip(lower=0))
        for feat in ['raw_ani', 'mash_ani', 'shared_log', 'af_q', 'af_r']:
            col = f'{prefix}_{feat}'
            df[col] = df[col].fillna(df[col].median())
            if feat in ('af_q', 'af_r'):
                df[col] = df[col].clip(0, 1)
            feature_cols.append(col)

    X_df = df[feature_cols].fillna(df[feature_cols].median())
    return X_df.values, df, base_prefix, feature_cols


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', required=True, help='Path to GBRT v5 pickle model')
    parser.add_argument('--matrix', required=True, help='Validation matrix TSV')
    parser.add_argument('--mode', required=True, choices=['bcgi', 'cjepi', 'combined'])
    parser.add_argument('--output', required=True, help='Output TSV path')
    parser.add_argument('--label', default='validation')
    args = parser.parse_args()

    with open(args.model, 'rb') as f:
        model = pickle.load(f)

    df = pd.read_csv(args.matrix, sep='\t', low_memory=False)
    X, df, base_prefix, feature_cols = build_features(df, args.mode)

    # Verify model expects the same number of features
    if hasattr(model, 'n_features_in_') and model.n_features_in_ != X.shape[1]:
        raise ValueError(
            f"Model expects {model.n_features_in_} features but got {X.shape[1]}. "
            f"Expected: {feature_cols}"
        )

    df['gbrt_v5_bias'] = model.predict(X)
    df['gbrt_v5_ani'] = (df['base_raw_ani'] + df['gbrt_v5_bias']).clip(0, 1)

    ref_col = 'fastani_ani'
    if ref_col in df.columns:
        mask = df[ref_col].notna() & df['gbrt_v5_ani'].notna()
        sub = df[mask].copy()

        methods = {
            f'Syn2bANI {base_prefix} raw': f'{base_prefix}_raw_ani',
            f'Syn2bANI {base_prefix} mash': f'{base_prefix}_mash_ani',
            'Syn2bANI GBRT v5': 'gbrt_v5_ani',
        }
        if 'skani_ani' in sub.columns:
            methods['skani'] = 'skani_ani'

        print(f"\n=== {args.label} Metrics ({args.mode}, n={len(sub)}) ===")
        print(f"{'Method':<30} {'MAE':>10} {'RMSE':>10}")
        rows = []
        for name, col in methods.items():
            if col in sub.columns:
                m = sub[col]
                r = sub[ref_col]
                mae = mean_absolute_error(r, m) * 100
                rmse = np.sqrt(mean_squared_error(r, m)) * 100
                print(f"{name:<30} {mae:>10.3f}% {rmse:>10.3f}%")
                rows.append({'method': name, 'mae': mae, 'rmse': rmse})

    df.to_csv(args.output, sep='\t', index=False, float_format='%.6f')
    print(f"\nSaved: {args.output}")


if __name__ == '__main__':
    main()
