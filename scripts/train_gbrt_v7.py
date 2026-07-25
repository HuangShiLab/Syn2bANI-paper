#!/usr/bin/env python3
"""Train GBRT v7 debiasing model with chain k-mer ANI feature.

Features: raw_ani, mash_ani, chained_kmer_ani, shared_log, af_q, af_r
Target:   fastani_ani - s2b_raw_ani
"""
import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error


def compute_mash_ani(af_q, af_r, tag_len=32.0):
    """Compute Mash-like ANI from bidirectional containment."""
    c1 = np.maximum(af_q, 1e-10)
    c2 = np.maximum(af_r, 1e-10)
    containment_geo = np.sqrt(c1 * c2)
    return np.clip(1.0 + np.log(containment_geo) / tag_len, 0.0, 1.0)


def build_features(df: pd.DataFrame, reference_col: str = 'fastani_ani'):
    """Build feature matrix X and target y for GBRT v7."""
    if reference_col not in df.columns:
        raise ValueError(f"Reference column '{reference_col}' not found.")

    # Keep rows with valid reference and raw ANI
    df = df[df[reference_col].notna() & df['s2b_raw_ani'].notna()].copy()

    # Compute mash_ani if not present
    if 's2b_mash_ani' not in df.columns:
        df['s2b_mash_ani'] = compute_mash_ani(df['s2b_af_q'].values, df['s2b_af_r'].values)

    # Target: bias to add to raw_ani
    df['y'] = df[reference_col] - df['s2b_raw_ani']

    feature_cols = []

    # 1. Raw ANI
    df['raw_ani'] = df['s2b_raw_ani']
    feature_cols.append('raw_ani')

    # 2. Mash-like ANI
    df['mash_ani'] = df['s2b_mash_ani']
    feature_cols.append('mash_ani')

    # 3. Chain-interval k-mer ANI
    df['chained_kmer_ani'] = df['s2b_chained_kmer_ani'].fillna(df['s2b_mash_ani']).clip(0, 1)
    feature_cols.append('chained_kmer_ani')

    # 4. Shared tag count (log1p)
    df['shared_log'] = np.log1p(df['s2b_shared_tags'].fillna(0).clip(lower=0))
    feature_cols.append('shared_log')

    # 5. Alignment fractions
    for suffix in ['_q', '_r']:
        col = f's2b_af{suffix}'
        df[f'af{suffix}'] = df[col].fillna(df[col].median()).clip(lower=0.0, upper=1.0)
        feature_cols.append(f'af{suffix}')

    X_df = df[feature_cols].copy()
    X_df = X_df.fillna(X_df.median())
    X = X_df.values
    y = df['y'].values

    level_col = 'label' if 'label' in df.columns else 'level'
    meta_cols = ['query', 'reference', reference_col, 's2b_raw_ani', 's2b_mash_ani', 's2b_chained_kmer_ani']
    if level_col in df.columns:
        meta_cols.append(level_col)
    meta = df[[c for c in meta_cols if c in df.columns]].copy()

    return X, y, meta, feature_cols, level_col


def train_model(X, y, feature_names):
    """Train GBRT with hyperparameter grid search."""
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )

    configs = [
        {'n_estimators': 300, 'max_depth': 5, 'learning_rate': 0.05, 'subsample': 0.8},
        {'n_estimators': 500, 'max_depth': 6, 'learning_rate': 0.03, 'subsample': 0.8},
        {'n_estimators': 500, 'max_depth': 5, 'learning_rate': 0.05, 'subsample': 0.9},
        {'n_estimators': 1000, 'max_depth': 5, 'learning_rate': 0.02, 'subsample': 0.8},
        {'n_estimators': 1000, 'max_depth': 6, 'learning_rate': 0.02, 'subsample': 0.8},
    ]

    best_mae = float('inf')
    best_model = None
    best_params = {}

    print("\nTraining configurations:")
    for cfg in configs:
        model = GradientBoostingRegressor(
            n_estimators=cfg['n_estimators'],
            max_depth=cfg['max_depth'],
            learning_rate=cfg['learning_rate'],
            subsample=cfg['subsample'],
            random_state=42,
            loss='squared_error'
        )
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        mae = mean_absolute_error(y_test, y_pred)
        print(f"  {cfg} -> Test MAE (bias): {mae:.4f}")

        if mae < best_mae:
            best_mae = mae
            best_model = model
            best_params = cfg

    print(f"\nBest config: {best_params} (bias MAE={best_mae:.4f})")

    # 5-fold CV
    print("\n5-fold cross-validation:")
    cv_scores = cross_val_score(
        best_model, X, y, cv=5, scoring='neg_mean_absolute_error'
    )
    print(f"  CV MAE (bias): {-cv_scores.mean():.4f} (+/- {cv_scores.std():.4f})")

    # Final fit on all data
    best_model.fit(X, y)

    print("\nFeature importances:")
    for name, imp in zip(feature_names, best_model.feature_importances_):
        print(f"  {name}: {imp:.4f}")

    return best_model, best_params


def _sklearn_tree_to_json(tree, feature_names):
    """Convert sklearn Tree to Rust GbrtModel format."""
    n_nodes = tree.node_count
    nodes = []
    for i in range(n_nodes):
        left = int(tree.children_left[i])
        right = int(tree.children_right[i])
        if left == right:
            value = float(tree.value[i][0][0])
            nodes.append({'type': 'leaf', 'value': value})
        else:
            feature = int(tree.feature[i])
            threshold = float(tree.threshold[i])
            nodes.append({
                'type': 'split',
                'feature': feature,
                'feature_name': feature_names[feature],
                'threshold': threshold,
                'left': left,
                'right': right,
            })
    return {'nodes': nodes}


def export_model(model, feature_names, output_path):
    """Export sklearn GBRT to Rust JSON format."""
    init_value = float(model.init_.constant_[0][0]) if hasattr(model, 'init_') else 0.0

    trees = []
    for est in model.estimators_.flatten():
        trees.append(_sklearn_tree_to_json(est.tree_, feature_names))

    model_data = {
        'meta': {
            'n_estimators': int(model.n_estimators),
            'max_depth': int(model.max_depth),
            'learning_rate': float(model.learning_rate),
            'init_value': init_value,
            'feature_names': list(feature_names),
        },
        'trees': trees,
    }

    with open(output_path, 'w') as f:
        json.dump(model_data, f, indent=2)

    pickle_path = str(output_path).replace('.json', '.pkl')
    with open(pickle_path, 'wb') as f:
        pickle.dump(model, f)

    print(f"\nModel exported:\n  JSON: {output_path}\n  Pickle: {pickle_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix', required=True)
    parser.add_argument('--output', default='results/gbrt_model_v7.json')
    parser.add_argument('--report', default='results/gbrt_v7_report.txt')
    args = parser.parse_args()

    print(f"Loading matrix: {args.matrix}")
    df = pd.read_csv(args.matrix, sep='\t', low_memory=False)
    print(f"Total rows: {len(df)}")

    X, y, meta, feature_names, level_col = build_features(df)
    print(f"Training samples: {len(y)}")
    print(f"Features: {feature_names}")

    model, params = train_model(X, y, feature_names)

    # In-sample predictions
    meta['predicted_bias'] = model.predict(X)
    meta['gbrt_v7_ani'] = (meta['s2b_raw_ani'] + meta['predicted_bias']).clip(0, 1)

    ref_col = 'fastani_ani'
    mask = meta[ref_col].notna() & meta['gbrt_v7_ani'].notna()
    sub = meta[mask]
    mae = mean_absolute_error(sub[ref_col], sub['gbrt_v7_ani'])
    rmse = np.sqrt(mean_squared_error(sub[ref_col], sub['gbrt_v7_ani']))
    print(f"\nIn-sample corrected ANI MAE: {mae*100:.4f}%")
    print(f"In-sample corrected ANI RMSE: {rmse*100:.4f}%")

    # Export
    export_model(model, feature_names, args.output)

    # Report
    lines = [
        '# GBRT v7 Training Report',
        '',
        f'Matrix: {args.matrix}',
        f'Training samples: {len(y)}',
        f'Features: {feature_names}',
        f'Best params: {params}',
        f'In-sample corrected ANI MAE: {mae*100:.4f}%',
        f'In-sample corrected ANI RMSE: {rmse*100:.4f}%',
        '',
        '## Feature importances',
    ]
    for name, imp in zip(feature_names, model.feature_importances_):
        lines.append(f'- {name}: {imp:.4f}')

    Path(args.report).write_text('\n'.join(lines))
    print(f"Report saved: {args.report}")


if __name__ == '__main__':
    main()
