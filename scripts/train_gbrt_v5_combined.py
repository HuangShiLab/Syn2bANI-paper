#!/usr/bin/env python3
"""Train GBRT v5 debiasing model with BcgI, CjePI, or combined features."""
import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error


def build_features(df, mode):
    """Build feature matrix X and target y."""
    df = df[df['fastani_ani'].notna()].copy()

    prefixes = []
    if mode in ('bcgi', 'combined'):
        prefixes.append('bcgi')
    if mode in ('cjepi', 'combined'):
        prefixes.append('cjepi')

    # Base estimator is raw_ani from the first available enzyme
    base_prefix = prefixes[0]
    df['base_raw_ani'] = df[f'{base_prefix}_raw_ani']

    # Target: bias to add to base raw_ani
    df['y'] = df['fastani_ani'] - df['base_raw_ani']

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
    X = X_df.values
    y = df['y'].values

    meta = df[['query', 'reference', 'fastani_ani', 'base_raw_ani']].copy()
    if 'label' in df.columns:
        meta['label'] = df['label']

    return X, y, meta, feature_cols


def train_model(X, y, feature_names):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

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
        print(f"  {cfg} -> Test bias MAE: {mae:.4f}")

        if mae < best_mae:
            best_mae = mae
            best_model = model
            best_params = cfg

    print(f"\nBest config: {best_params} (bias MAE={best_mae:.4f})")

    print("\n5-fold cross-validation:")
    cv_scores = cross_val_score(best_model, X, y, cv=5, scoring='neg_mean_absolute_error')
    print(f"  CV bias MAE: {-cv_scores.mean():.4f} (+/- {cv_scores.std():.4f})")

    best_model.fit(X, y)

    print("\nFeature importances:")
    for name, imp in zip(feature_names, best_model.feature_importances_):
        print(f"  {name}: {imp:.4f}")

    return best_model, best_params


def _sklearn_tree_to_json(tree, feature_names):
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
    parser.add_argument('--mode', required=True, choices=['bcgi', 'cjepi', 'combined'])
    parser.add_argument('--output', required=True)
    parser.add_argument('--report', required=True)
    args = parser.parse_args()

    print(f"Loading matrix: {args.matrix}")
    df = pd.read_csv(args.matrix, sep='\t', low_memory=False)
    print(f"Total rows: {len(df)}")

    X, y, meta, feature_names = build_features(df, args.mode)
    print(f"Training samples: {len(y)}")
    print(f"Features: {feature_names}")

    model, params = train_model(X, y, feature_names)

    meta['predicted_bias'] = model.predict(X)
    meta['gbrt_v5_ani'] = (meta['base_raw_ani'] + meta['predicted_bias']).clip(0, 1)

    mask = meta['fastani_ani'].notna() & meta['gbrt_v5_ani'].notna()
    sub = meta[mask]
    mae = mean_absolute_error(sub['fastani_ani'], sub['gbrt_v5_ani'])
    rmse = np.sqrt(mean_squared_error(sub['fastani_ani'], sub['gbrt_v5_ani']))
    print(f"\nIn-sample corrected ANI MAE: {mae*100:.4f}%")
    print(f"In-sample corrected ANI RMSE: {rmse*100:.4f}%")

    export_model(model, feature_names, args.output)

    lines = [
        f'# GBRT v5 {args.mode} Training Report',
        '',
        f'Matrix: {args.matrix}',
        f'Mode: {args.mode}',
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
