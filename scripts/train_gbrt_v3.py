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


def build_features(df: pd.DataFrame, use_skani_features: bool = False,
                   reference_col: str = 'fastani_ani') -> tuple:
    """Build feature matrix X and target y from benchmark results.

    By default only features available at Syn2bANI inference time are used
    (raw ANI, shared tag count, taxonomy label, GC, genome size). skani-derived
    features are opt-in because they are not available when Syn2bANI runs alone.
    """
    if reference_col not in df.columns:
        raise ValueError(f"Reference column '{reference_col}' not found in matrix. "
                         f"Available: {list(df.columns)}")

    # Filter rows with valid reference ANI and raw Syn2bANI output
    df = df[df[reference_col].notna() & df['s2b_raw_ani'].notna()].copy()

    # Target: the bias that GBRT should learn to correct
    # Positive = Syn2bANI underestimates, Negative = overestimates
    df['y'] = df[reference_col] - df['s2b_raw_ani']

    # Features (names must match the Rust inference code)
    feature_cols = []

    # 1. Raw ANI
    df['raw_ani'] = df['s2b_raw_ani']
    feature_cols.append('raw_ani')

    # 2. Shared tag count (log1p transformed)
    if 's2b_shared_tags' in df.columns and df['s2b_shared_tags'].notna().any():
        df['shared_log'] = np.log1p(df['s2b_shared_tags'].fillna(0).clip(lower=0))
        feature_cols.append('shared_log')

    # 3. Taxonomy relatedness label (used for reporting; optional as a model feature
    # because it is not available at inference time unless taxonomy is supplied).
    use_taxonomy_feature = False  # set to True to include label as a feature
    if 'level' in df.columns:
        level_col = 'level'
        level_map = {'intra_species': 0, 'intra_genus': 1,
                     'intra_family': 2, 'random': 3}
    else:
        level_col = 'label'
        level_map = {'high': 0, 'mid_high': 1, 'mid': 2, 'low': 3}
    if use_taxonomy_feature:
        df['level'] = df[level_col].map(level_map).fillna(3)
        feature_cols.append('level')

    # 4. GC content (if available from metadata)
    for suffix in ['_q', '_r']:
        col = f'gc_percentage{suffix}'
        if col in df.columns:
            df[f'gc{suffix}'] = df[col].fillna(df[col].median())
            feature_cols.append(f'gc{suffix}')

    # 5. Syn2bANI alignment fractions (inference-time features)
    for suffix in ['_q', '_r']:
        col = f's2b_af{suffix}'
        if col in df.columns:
            df[f'af{suffix}'] = df[col].fillna(df[col].median()).clip(lower=0.0, upper=1.0)
            feature_cols.append(f'af{suffix}')

    # 6. Reference GC content (inference-time feature, but not passed to
    # AniCalculator in the normal dist flow, so exclude from production model).
    # if 's2b_ref_gc' in df.columns:
    #     df['ref_gc'] = df['s2b_ref_gc'].fillna(df['s2b_ref_gc'].median())
    #     feature_cols.append('ref_gc')

    # 7. Genome size (if available)
    for suffix in ['_q', '_r']:
        col = f'genome_size{suffix}'
        if col in df.columns:
            df[f'genome_size{suffix}'] = np.log10(df[col].fillna(df[col].median()) + 1)
            feature_cols.append(f'genome_size{suffix}')

    # 8. skani alignment fraction (only if explicitly requested — leaks skani info)
    if use_skani_features and 'skani_align_frac' in df.columns:
        df['skani_align_frac'] = df['skani_align_frac'].fillna(0.5)
        feature_cols.append('skani_align_frac')

    X_df = df[feature_cols].copy()
    X_df = X_df.fillna(X_df.median())
    X = X_df.values
    y = df['y'].values
    meta = df[['query', 'reference', level_col,
               reference_col, 'skani_ani',
               's2b_raw_ani', 's2b_gbrt_ani', 'q_phylum', 'r_phylum']].copy()

    return X, y, meta, feature_cols, level_col


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


def _sklearn_tree_to_json(tree, feature_names):
    """Convert one sklearn Tree object to the Rust GbrtModel tree format."""
    n_nodes = tree.node_count
    nodes = []
    for i in range(n_nodes):
        left = int(tree.children_left[i])
        right = int(tree.children_right[i])
        if left == right:  # leaf
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
    """Export sklearn GBRT to the Rust GbrtModel JSON format."""
    # Initial value for squared-error loss is the mean target.
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

    # Also save pickle for Python reuse / inspection
    pickle_path = str(output_path).replace('.json', '.pkl')
    with open(pickle_path, 'wb') as f:
        pickle.dump(model, f)

    print(f"\nModel exported:")
    print(f"  Rust JSON: {output_path}")
    print(f"  Pickle:    {pickle_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix', required=True, help='Benchmark matrix TSV')
    parser.add_argument('--output', default='results/gbrt_model_v3.json')
    parser.add_argument('--report', default='results/gbrt_v3_report.txt')
    parser.add_argument('--reference', default='fastani_ani',
                        choices=['fastani_ani', 'skani_ani'],
                        help='Reference ANI column to train against (default: fastani_ani)')
    parser.add_argument('--use-skani-features', action='store_true',
                        help='Allow skani_align_frac as a feature (leaks skani info)')
    args = parser.parse_args()

    print(f"Loading matrix: {args.matrix}")
    df = pd.read_csv(args.matrix, sep='\t', low_memory=False)
    print(f"Total rows: {len(df)}")

    X, y, meta, feature_names, level_col = build_features(
        df,
        use_skani_features=args.use_skani_features,
        reference_col=args.reference,
    )
    print(f"Training samples: {len(y)}")
    print(f"Features: {feature_names}")
    print(f"Reference: {args.reference}")
    print(f"Target range: [{y.min():.4f}, {y.max():.4f}]")
    print(f"Target mean: {y.mean():.4f}, std: {y.std():.4f}")

    model, params = train_model(X, y, feature_names)

    # Predict on all data
    y_pred = model.predict(X)
    meta['predicted_bias'] = y_pred
    meta['corrected_ani'] = meta['s2b_raw_ani'] + y_pred

    # Evaluate against the same reference used for training
    ref = meta[args.reference]
    mae_before = mean_absolute_error(ref, meta['s2b_raw_ani'])
    mae_after = mean_absolute_error(ref, meta['corrected_ani'])
    rmse_before = np.sqrt(mean_squared_error(ref, meta['s2b_raw_ani']))
    rmse_after = np.sqrt(mean_squared_error(ref, meta['corrected_ani']))

    ref_name = 'FastANI' if args.reference == 'fastani_ani' else 'skani'
    report = f"""
GBRT v3 Training Report
=======================
Training samples: {len(y)}
Features: {feature_names}
Reference: {args.reference}
Best params: {params}

MAE (raw Syn2bANI vs {ref_name}):     {mae_before:.4f}
MAE (GBRT corrected vs {ref_name}):   {mae_after:.4f}
Improvement:                          {mae_before - mae_after:.4f} ({(1 - mae_after/mae_before)*100:.1f}%)

RMSE (raw):                           {rmse_before:.4f}
RMSE (GBRT):                          {rmse_after:.4f}

By relatedness label:
"""
    for level in sorted(meta[level_col].unique()):
        subset = meta[meta[level_col] == level]
        mae_b = mean_absolute_error(subset[args.reference], subset['s2b_raw_ani'])
        mae_a = mean_absolute_error(subset[args.reference], subset['corrected_ani'])
        report += f"  {level:20s}: raw={mae_b:.4f}, gbrt={mae_a:.4f} (n={len(subset)})\n"

    # Per-phylum error report using query phylum as the grouping key
    meta['q_phylum_clean'] = meta['q_phylum'].str.replace(r'^p__', '', regex=True)
    report += f"\nPer-phylum error (GBRT corrected vs {ref_name}, top 15 by count):\n"
    phylum_counts = meta['q_phylum_clean'].value_counts()
    for phylum in phylum_counts.head(15).index:
        subset = meta[meta['q_phylum_clean'] == phylum]
        mae_b = mean_absolute_error(subset[args.reference], subset['s2b_raw_ani'])
        mae_a = mean_absolute_error(subset[args.reference], subset['corrected_ani'])
        report += f"  {phylum:30s}: raw={mae_b:.4f}, gbrt={mae_a:.4f} (n={len(subset)})\n"

    print(report)

    with open(args.report, 'w') as f:
        f.write(report)
    print(f"Report saved: {args.report}")

    export_model(model, feature_names, args.output)


if __name__ == '__main__':
    main()
