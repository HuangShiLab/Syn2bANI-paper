#!/usr/bin/env python3
"""Evaluate Syn2bANI GBRT v4 on independent oral/gut validation set."""
import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error, mean_squared_error

import sys
sys.path.insert(0, str(Path(__file__).parent))
from train_gbrt_v3 import build_features


def compute_metrics(df, pred_col, ref_col='fastani_ani'):
    mask = df[ref_col].notna() & df[pred_col].notna()
    sub = df[mask]
    if len(sub) == 0:
        return None
    mae = mean_absolute_error(sub[ref_col], sub[pred_col])
    rmse = np.sqrt(mean_squared_error(sub[ref_col], sub[pred_col]))
    r, _ = pearsonr(sub[ref_col], sub[pred_col])
    return {'n': len(sub), 'mae': mae, 'rmse': rmse, 'r': r}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix', required=True)
    parser.add_argument('--model', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--report', required=True)
    args = parser.parse_args()

    # Load model
    with open(args.model, 'rb') as f:
        model = pickle.load(f)

    # Load validation matrix
    df = pd.read_csv(args.matrix, sep='\t', low_memory=False)
    print(f'Loaded validation matrix: {len(df)} rows')

    # Add dummy phylum columns expected by build_features
    for col in ['q_phylum', 'r_phylum']:
        if col not in df.columns:
            df[col] = 'unknown'

    # Build features same as training
    X, y, meta, feature_names, level_col = build_features(
        df, use_skani_features=False, reference_col='fastani_ani'
    )
    print(f'Features: {feature_names}')
    print(f'Rows with FastANI reference: {len(meta)}')

    # Predict GBRT v4 bias and corrected ANI
    predicted_bias = model.predict(X)
    meta['gbrt_v4_bias'] = predicted_bias
    meta['gbrt_v4_ani'] = meta['s2b_raw_ani'] + predicted_bias

    # Merge back into df
    df = df.merge(
        meta[['query', 'reference', 'gbrt_v4_ani']],
        on=['query', 'reference'], how='left'
    )

    # Compute metrics
    metrics = {}
    metrics['Syn2bANI raw'] = compute_metrics(df, 's2b_raw_ani')
    metrics['Syn2bANI GBRT v4'] = compute_metrics(df, 'gbrt_v4_ani')
    metrics['skani'] = compute_metrics(df, 'skani_ani')
    metrics['Syn2bANI built-in GBRT'] = compute_metrics(df, 's2b_gbrt_ani')

    # Per-label metrics
    per_label = {}
    for label in df[level_col].unique():
        if pd.isna(label):
            continue
        sub = df[df[level_col] == label]
        per_label[label] = {
            'Syn2bANI raw': compute_metrics(sub, 's2b_raw_ani'),
            'Syn2bANI GBRT v4': compute_metrics(sub, 'gbrt_v4_ani'),
            'skani': compute_metrics(sub, 'skani_ani'),
            'Syn2bANI built-in GBRT': compute_metrics(sub, 's2b_gbrt_ani'),
        }

    # Write report
    lines = []
    lines.append('# Independent Oral/Gut Validation Report')
    lines.append('')
    lines.append(f'Validation matrix: `{args.matrix}`')
    lines.append(f'Trained model: `{args.model}`')
    lines.append(f'Total pairs: {len(df)}')
    lines.append(f'Pairs with FastANI reference: {df["fastani_ani"].notna().sum()}')
    lines.append('')

    lines.append('## Overall Metrics (vs FastANI reference)')
    lines.append('')
    lines.append('| Method | n | MAE | RMSE | Pearson r |')
    lines.append('|--------|---|-----|------|-----------|')
    for method, m in metrics.items():
        if m is None:
            lines.append(f'| {method} | - | - | - | - |')
        else:
            lines.append(f'| {method} | {m["n"]} | {m["mae"]*100:.3f}% | {m["rmse"]*100:.3f}% | {m["r"]:.4f} |')
    lines.append('')

    lines.append('## Per-Label Metrics')
    lines.append('')
    for label, methods in per_label.items():
        lines.append(f'### {label}')
        lines.append('')
        lines.append('| Method | n | MAE | RMSE | Pearson r |')
        lines.append('|--------|---|-----|------|-----------|')
        for method, m in methods.items():
            if m is None:
                lines.append(f'| {method} | - | - | - | - |')
            else:
                lines.append(f'| {method} | {m["n"]} | {m["mae"]*100:.3f}% | {m["rmse"]*100:.3f}% | {m["r"]:.4f} |')
        lines.append('')

    Path(args.report).write_text('\n'.join(lines))
    print(f'Report saved: {args.report}')

    # Save updated matrix with GBRT v4 predictions
    df.to_csv(args.output, sep='\t', index=False, float_format='%.6f')
    print(f'Updated matrix saved: {args.output}')

    # Print summary
    print('\n=== Overall Metrics ===')
    for method, m in metrics.items():
        if m:
            print(f'{method}: n={m["n"]}, MAE={m["mae"]*100:.3f}%, RMSE={m["rmse"]*100:.3f}%, r={m["r"]:.4f}')


if __name__ == '__main__':
    main()
