#!/usr/bin/env python3
"""
train_gbrt_v3.py
Train GBRT v3 debiasing model from benchmark matrix results.

Usage:
  python3 train_gbrt_v3.py \
    --matrix results/matrix.tsv \
    --output results/gbrt_model_v3.json \
    --report results/gbrt_v3_report.txt
"""

import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


def build_features(df: pd.DataFrame) -> tuple:
    """Build feature matrix X and target y from benchmark results."""
    # Filter rows with valid reference ANI (skani as ground truth)
    df = df[df['skani_ani'].notna() & df['s2b_raw_ani'].notna()].copy()

    # Target: the bias that GBRT should learn to correct
    # Positive = Syn2bANI underestimates, Negative = overestimates
    df['y'] = df['skani_ani'] - df['s2b_raw_ani']

    # Features
    feature_cols = []

    # 1. Raw ANI features
    df['feat_ani_raw'] = df['s2b_raw_ani']
    feature_cols.append('feat_ani_raw')

    # 2. Shared tag ratio (if available)
    if 's2b_shared_tags' in df.columns and df['s2b_shared_tags'].notna().any():
        # Normalize by some proxy for total tags (use median if unavailable)
        df['feat_shared_ratio'] = df['s2b_shared_tags'] / df['s2b_shared_tags'].median()
        feature_cols.append('feat_shared_ratio')

    # 3. Taxonomy level (categorical → numeric)
    level_map = {'intra_species': 0, 'intra_genus': 1,
                 'intra_family': 2, 'random': 3}
    df['feat_level'] = df['level'].map(level_map).fillna(3)
    feature_cols.append('feat_level')

    # 4. GC content (if available from metadata)
    for suffix in ['_q', '_r']:
        col = f'gc_percentage{suffix}'
        if col in df.columns:
            df[f'feat_gc{suffix}'] = df[col].fillna(df[col].median())
            feature_cols.append(f'feat_gc{suffix}')

    # 5. Genome size (if available)
    for suffix in ['_q', '_r']:
        col = f'genome_size{suffix}'
        if col in df.columns:
            df[f'feat_size{suffix}'] = np.log10(df[col].fillna(df[col].median()) + 1)
            feature_cols.append(f'feat_size{suffix}')

    # 6. skani alignment fraction (if available)
    if 'skani_align_frac' in df.columns:
        df['feat_align_frac'] = df['skani_align_frac'].fillna(0.5)
        feature_cols.append('feat_align_frac')

    X = df[feature_cols].values
    y = df['y'].values
    meta = df[['query', 'reference', 'level', 'skani_ani',
               's2b_raw_ani', 's2b_gbrt_ani']].copy()

    return X, y, meta, feature_cols


def train_model(X, y, feature_names):
    """Train GBRT model with cross-validation."""
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=None
    )

    # Grid search over key hyperparameters
    best_mae = float('inf')
    best_model = None
    best_params = {}

    configs = [
        {'n_estimators': 300, 'max_depth': 5, 'learning_rate': 0.05, 'subsample': 0.8},
        {'n_estimators': 500, 'max_depth': 6, 'learning_rate': 0.03, 'subsample': 0.8},
        {'n_estimators': 500, 'max_depth': 5, 'learning_rate': 0.05, 'subsample': 0.9},
        {'n_estimators': 1000, 'max_depth': 5, 'learning_rate': 0.02, 'subsample': 0.8},
    ]

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
        print(f"  {cfg} -> Test MAE: {mae:.4f}")

        if mae < best_mae:
            best_mae = mae
            best_model = model
            best_params = cfg

    print(f"\nBest config: {best_params} (MAE={best_mae:.4f})")

    # Cross-validation on full data
    print("\n5-fold cross-validation:")
    cv_scores = cross_val_score(
        best_model, X, y, cv=5,
        scoring='neg_mean_absolute_error'
    )
    print(f"  CV MAE: {-cv_scores.mean():.4f} (+/- {cv_scores.std():.4f})")

    # Final fit on all data
    best_model.fit(X, y)

    # Feature importance
    print("\nFeature importances:")
    for name, imp in zip(feature_names, best_model.feature_importances_):
        print(f"  {name}: {imp:.4f}")

    return best_model, best_params


def export_model(model, feature_names, output_path):
    """Export model to JSON format for embedding in Rust binary."""
    # sklearn GBRT doesn't have a simple JSON export, so we use pickle
    # For Rust embedding, we need a custom serialization

    # Simplified: store key parameters + trees as nested dicts
    # Full tree serialization is complex; for production, consider:
    #  - treelite (https://treelite.readthedocs.io/)
    #  - ONNX export
    #  - Custom JSON tree walker

    model_data = {
        'version': '3.0',
        'n_estimators': model.n_estimators,
        'max_depth': model.max_depth,
        'learning_rate': model.learning_rate,
        'subsample': model.subsample,
        'feature_names': feature_names,
        'feature_importances': model.feature_importances_.tolist(),
        'train_score': model.train_score_.tolist() if hasattr(model, 'train_score_') else [],
        'n_features': model.n_features_in_,
        # Tree structure is too complex for simple JSON;
        # use pickle for now, convert to custom format later
    }

    with open(output_path, 'w') as f:
        json.dump(model_data, f, indent=2)

    # Also save pickle for Python reuse
    pickle_path = str(output_path).replace('.json', '.pkl')
    with open(pickle_path, 'wb') as f:
        pickle.dump(model, f)

    print(f"\nModel exported:")
    print(f"  JSON (params): {output_path}")
    print(f"  Pickle (full): {pickle_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix', required=True, help='Benchmark matrix TSV')
    parser.add_argument('--output', default='results/gbrt_model_v3.json')
    parser.add_argument('--report', default='results/gbrt_v3_report.txt')
    args = parser.parse_args()

    print(f"Loading matrix: {args.matrix}")
    df = pd.read_csv(args.matrix, sep='\t', low_memory=False)
    print(f"Total rows: {len(df)}")

    X, y, meta, feature_names = build_features(df)
    print(f"Training samples: {len(y)}")
    print(f"Features: {feature_names}")
    print(f"Target range: [{y.min():.4f}, {y.max():.4f}]")
    print(f"Target mean: {y.mean():.4f}, std: {y.std():.4f}")

    model, params = train_model(X, y, feature_names)

    # Predict on all data
    y_pred = model.predict(X)
    meta['predicted_bias'] = y_pred
    meta['corrected_ani'] = meta['s2b_raw_ani'] + y_pred

    # Evaluate
    mae_before = mean_absolute_error(meta['skani_ani'], meta['s2b_raw_ani'])
    mae_after = mean_absolute_error(meta['skani_ani'], meta['corrected_ani'])
    rmse_before = np.sqrt(mean_squared_error(meta['skani_ani'], meta['s2b_raw_ani']))
    rmse_after = np.sqrt(mean_squared_error(meta['skani_ani'], meta['corrected_ani']))

    report = f"""
GBRT v3 Training Report
=======================
Training samples: {len(y)}
Features: {feature_names}
Best params: {params}

MAE (raw Syn2bANI vs skani):     {mae_before:.4f}
MAE (GBRT corrected vs skani):   {mae_after:.4f}
Improvement:                     {mae_before - mae_after:.4f} ({(1 - mae_after/mae_before)*100:.1f}%)

RMSE (raw):                      {rmse_before:.4f}
RMSE (GBRT):                     {rmse_after:.4f}

By taxonomic level:
"""
    for level in meta['level'].unique():
        subset = meta[meta['level'] == level]
        mae_b = mean_absolute_error(subset['skani_ani'], subset['s2b_raw_ani'])
        mae_a = mean_absolute_error(subset['skani_ani'], subset['corrected_ani'])
        report += f"  {level:20s}: raw={mae_b:.4f}, gbrt={mae_a:.4f}\n"

    print(report)

    with open(args.report, 'w') as f:
        f.write(report)
    print(f"Report saved: {args.report}")

    export_model(model, feature_names, args.output)


if __name__ == '__main__':
    main()
