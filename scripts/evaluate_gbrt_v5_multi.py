#!/usr/bin/env python3
"""Evaluate multi-enzyme GBRT v5 model on a validation matrix."""
import argparse
import pickle

import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error, mean_squared_error


def build_features(df):
    df = df[df['fastani_ani'].notna()].copy()
    df['base_raw_ani'] = df['multi_raw_ani']
    df['multi_shared_log'] = np.log1p(df['multi_shared_tags'].fillna(0).clip(lower=0))
    feature_cols = ['multi_raw_ani', 'multi_mash_ani', 'multi_shared_log', 'multi_af_q', 'multi_af_r']
    for col in feature_cols:
        df[col] = df[col].fillna(df[col].median())
        if col in ('multi_af_q', 'multi_af_r'):
            df[col] = df[col].clip(0, 1)
    X_df = df[feature_cols].fillna(df[feature_cols].median())
    return X_df.values, df, feature_cols


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
    X, df, feature_cols = build_features(df)

    if hasattr(model, 'n_features_in_') and model.n_features_in_ != X.shape[1]:
        raise ValueError(f"Model expects {model.n_features_in_} features but got {X.shape[1]}")

    df['gbrt_v5_multi_bias'] = model.predict(X)
    df['gbrt_v5_multi_ani'] = (df['base_raw_ani'] + df['gbrt_v5_multi_bias']).clip(0, 1)

    ref_col = 'fastani_ani'
    if ref_col in df.columns:
        mask = df[ref_col].notna() & df['gbrt_v5_multi_ani'].notna()
        sub = df[mask].copy()
        methods = {
            'Syn2bANI multi raw': 'multi_raw_ani',
            'Syn2bANI multi mash': 'multi_mash_ani',
            'Syn2bANI GBRT v5 multi': 'gbrt_v5_multi_ani',
        }
        if 'skani_ani' in sub.columns:
            methods['skani'] = 'skani_ani'

        print(f"\n=== {args.label} Metrics (multi, n={len(sub)}) ===")
        print(f"{'Method':<30} {'MAE':>10} {'RMSE':>10}")
        for name, col in methods.items():
            if col in sub.columns:
                m = sub[col]
                r = sub[ref_col]
                mae = mean_absolute_error(r, m) * 100
                rmse = np.sqrt(mean_squared_error(r, m)) * 100
                print(f"{name:<30} {mae:>10.3f}% {rmse:>10.3f}%")

    df.to_csv(args.output, sep='\t', index=False, float_format='%.6f')
    print(f"\nSaved: {args.output}")


if __name__ == '__main__':
    main()
