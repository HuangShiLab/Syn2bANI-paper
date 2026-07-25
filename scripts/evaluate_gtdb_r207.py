#!/usr/bin/env python3
"""
evaluate_gtdb_r207.py
Evaluate Syn2bANI (raw + GBRT) and skani against FastANI reference on GTDB-R207 pairs.

Generates:
  - results/gtdb_r207_evaluation.json
  - figures/gtdb_r207_scatter.png
  - figures/gtdb_r207_error_by_label.png
  - figures/gtdb_r207_phylum_error.png

Usage:
  python3 scripts/evaluate_gtdb_r207.py \
    --matrix results/matrix.tsv \
    --output results/gtdb_r207_evaluation.json \
    --figures figures/
"""

import argparse
import json
import os
import pickle
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def add_v4_predictions(df: pd.DataFrame, model_path: Path = None) -> pd.DataFrame:
    """Add s2b_v4_ani column using the clean GBRT v4 model pickle."""
    if model_path is None:
        root = Path(__file__).parent.parent
        model_path = root / 'results' / 'gbrt_model_v4_10k.pkl'
        if not model_path.exists():
            model_path = root / 'results' / 'gbrt_model_v4_1k.pkl'
    if not model_path.exists():
        return df
    try:
        with open(model_path, 'rb') as f:
            model = pickle.load(f)
    except Exception as e:
        print(f'WARNING: could not load v4 model {model_path}: {e}', file=sys.stderr)
        return df
    X = np.column_stack([
        df['s2b_raw_ani'].fillna(df['s2b_raw_ani'].median()).values,
        np.log1p(df['s2b_shared_tags'].fillna(0).clip(lower=0).values),
        df['s2b_af_q'].fillna(df['s2b_af_q'].median()).values,
        df['s2b_af_r'].fillna(df['s2b_af_r'].median()).values,
    ])
    bias = model.predict(X)
    df = df.copy()
    df['s2b_v4_ani'] = (df['s2b_raw_ani'].fillna(0) + bias).clip(0, 1)
    return df


def add_v7_predictions(df: pd.DataFrame, model_path: Path = None) -> pd.DataFrame:
    """Add s2b_v7_ani column using the GBRT v7 model pickle."""
    if model_path is None:
        root = Path(__file__).parent.parent
        model_path = root / 'results' / 'gbrt_model_v7_100k.pkl'
        if not model_path.exists():
            model_path = root / 'results' / 'gbrt_model_v7.pkl'
    if not model_path.exists():
        return df
    try:
        with open(model_path, 'rb') as f:
            model = pickle.load(f)
    except Exception as e:
        print(f'WARNING: could not load v7 model {model_path}: {e}', file=sys.stderr)
        return df
    raw_ani = df['s2b_raw_ani'].fillna(0)
    mash_ani = df['s2b_mash_ani'].combine_first(raw_ani)
    chained_kmer_ani = df['s2b_chained_kmer_ani'].combine_first(mash_ani)
    shared_log = np.log1p(df['s2b_shared_tags'].fillna(0).clip(lower=0).values)
    af_q_median = df['s2b_af_q'].median() if not pd.isna(df['s2b_af_q'].median()) else 0.0
    af_r_median = df['s2b_af_r'].median() if not pd.isna(df['s2b_af_r'].median()) else 0.0
    af_q = df['s2b_af_q'].fillna(af_q_median).values
    af_r = df['s2b_af_r'].fillna(af_r_median).values
    X = np.column_stack([raw_ani.values, mash_ani.values, chained_kmer_ani.values, shared_log, af_q, af_r])
    bias = model.predict(X)
    df = df.copy()
    df['s2b_v7_ani'] = (df['s2b_raw_ani'].fillna(0) + bias).clip(0, 1)
    return df


def compute_metrics(ref: np.ndarray, pred: np.ndarray, mask: np.ndarray = None) -> dict:
    """Compute MAE, RMSE, Pearson r."""
    if mask is not None:
        ref = ref[mask]
        pred = pred[mask]
    if len(ref) == 0:
        return {'n': 0, 'mae': None, 'rmse': None, 'pearson': None, 'mean_error': None}
    err = pred - ref
    mae = np.mean(np.abs(err))
    rmse = np.sqrt(np.mean(err ** 2))
    pearson = np.corrcoef(ref, pred)[0, 1] if len(ref) > 1 else 0.0
    return {
        'n': int(len(ref)),
        'mae': float(mae),
        'rmse': float(rmse),
        'pearson': float(pearson),
        'mean_error': float(np.mean(err)),
    }


def evaluate_by_group(df: pd.DataFrame, group_col: str, ref_col: str,
                      pred_cols: list[str]) -> dict:
    """Evaluate metrics per group."""
    results = {}
    groups = sorted(df[group_col].dropna().unique())
    for group in groups:
        mask = df[group_col] == group
        results[group] = {}
        for col in pred_cols:
            submask = mask & df[col].notna() & df[ref_col].notna()
            results[group][col] = compute_metrics(
                df[ref_col].values, df[col].values, submask.values
            )
    return results


def plot_scatter(df: pd.DataFrame, ref_col: str, pred_cols: dict[str, str],
                 out_path: Path):
    """Scatter plot of predictions vs reference."""
    n = len(pred_cols)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))
    if n == 1:
        axes = [axes]

    for ax, (name, col) in zip(axes, pred_cols.items()):
        sub = df.dropna(subset=[ref_col, col])
        x = sub[ref_col].values
        y = sub[col].values
        ax.scatter(x, y, alpha=0.5, s=20, edgecolors='white', linewidth=0.3)
        ax.plot([0.7, 1.0], [0.7, 1.0], 'k--', alpha=0.3, lw=1)
        ax.set_xlabel('FastANI (reference)', fontsize=11)
        ax.set_ylabel(name, fontsize=11)
        ax.set_xlim(0.7, 1.0)
        ax.set_ylim(0.7, 1.0)
        ax.grid(True, alpha=0.3)
        if len(x) > 0:
            mae = np.mean(np.abs(y - x))
            ax.set_title(f'{name}\nMAE = {mae*100:.2f}% (n={len(x)})', fontsize=12)

    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved scatter plot: {out_path}')


def plot_error_by_label(df: pd.DataFrame, ref_col: str, pred_cols: dict[str, str],
                        out_path: Path):
    """Boxplot of errors by label."""
    records = []
    for label in df['label'].unique():
        sub = df[df['label'] == label]
        for name, col in pred_cols.items():
            valid = sub.dropna(subset=[ref_col, col])
            for _, row in valid.iterrows():
                records.append({
                    'label': label,
                    'method': name,
                    'error': row[col] - row[ref_col],
                })

    if not records:
        print('No data for error-by-label plot')
        return

    err_df = pd.DataFrame(records)
    order = ['high', 'mid_high', 'mid', 'low']
    err_df['label'] = pd.Categorical(err_df['label'], categories=order, ordered=True)

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(data=err_df, x='label', y='error', hue='method', ax=ax)
    ax.axhline(0, color='red', linestyle='--', alpha=0.5)
    ax.set_xlabel('Pair label', fontsize=12)
    ax.set_ylabel('Error (predicted − FastANI)', fontsize=12)
    ax.set_title('ANI error distribution by pair type', fontsize=13)
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved error-by-label plot: {out_path}')


def plot_phylum_error(df: pd.DataFrame, ref_col: str, pred_cols: dict[str, str],
                      out_path: Path, min_pairs: int = 5):
    """Bar plot of MAE by query phylum."""
    phylum_stats = defaultdict(lambda: defaultdict(list))
    for _, row in df.iterrows():
        phylum = row.get('q_phylum', 'unknown')
        if pd.isna(phylum) or phylum == 'unknown':
            continue
        for name, col in pred_cols.items():
            if pd.notna(row[ref_col]) and pd.notna(row[col]):
                phylum_stats[phylum][name].append(abs(row[col] - row[ref_col]))

    phyla = []
    methods = list(pred_cols.keys())
    data = {m: [] for m in methods}
    for phylum in sorted(phylum_stats.keys()):
        n = min(len(phylum_stats[phylum][m]) for m in methods if phylum_stats[phylum][m])
        if n < min_pairs:
            continue
        phyla.append(phylum)
        for m in methods:
            vals = phylum_stats[phylum].get(m, [])
            data[m].append(np.mean(vals) if vals else np.nan)

    if not phyla:
        print('No phylum data for plot')
        return

    x = np.arange(len(phyla))
    width = 0.8 / len(methods)
    fig, ax = plt.subplots(figsize=(max(10, len(phyla) * 0.5), 6))
    for i, m in enumerate(methods):
        ax.bar(x + i * width - 0.4 + width / 2, data[m], width, label=m)

    ax.set_xticks(x)
    ax.set_xticklabels(phyla, rotation=45, ha='right')
    ax.set_ylabel('MAE vs FastANI', fontsize=12)
    ax.set_title(f'Per-phylum MAE (≥{min_pairs} pairs)', fontsize=13)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved phylum-error plot: {out_path}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--figures', required=True)
    parser.add_argument('--model', type=Path, default=None,
                        help='Path to GBRT v4 pickle (default: results/gbrt_model_v4_10k.pkl)')
    parser.add_argument('--model-v7', type=Path, default=None,
                        help='Path to GBRT v7 pickle (default: results/gbrt_model_v7_100k.pkl)')
    parser.add_argument('--min-phylum-pairs', type=int, default=5)
    args = parser.parse_args()

    df = pd.read_csv(args.matrix, sep='\t')
    print(f'Loaded {len(df)} rows from {args.matrix}')

    # Add v4 predictions if the model pickle and required columns are available
    if {'s2b_raw_ani', 's2b_shared_tags', 's2b_af_q', 's2b_af_r'}.issubset(df.columns):
        df = add_v4_predictions(df, model_path=args.model)

    # Add v7 predictions if chain k-mer features are available
    if {'s2b_raw_ani', 's2b_mash_ani', 's2b_chained_kmer_ani',
        's2b_shared_tags', 's2b_af_q', 's2b_af_r'}.issubset(df.columns):
        df = add_v7_predictions(df, model_path=args.model_v7)

    ref_col = 'fastani_ani'
    pred_cols = {
        'Syn2bANI raw': 's2b_raw_ani',
    }
    if 's2b_v7_ani' in df.columns:
        pred_cols['Syn2bANI GBRT v7'] = 's2b_v7_ani'
    elif 's2b_gbrt_ani' in df.columns:
        pred_cols['Syn2bANI GBRT'] = 's2b_gbrt_ani'
    if 's2b_v4_ani' in df.columns:
        pred_cols['Syn2bANI GBRT v4'] = 's2b_v4_ani'
    pred_cols['skani'] = 'skani_ani'

    # Filter to rows with reference
    df_valid = df[df[ref_col].notna()].copy()
    print(f'Rows with FastANI reference: {len(df_valid)}')

    # Global metrics
    global_metrics = {}
    for name, col in pred_cols.items():
        mask = df_valid[col].notna()
        global_metrics[name] = compute_metrics(
            df_valid[ref_col].values, df_valid[col].values, mask.values
        )

    # Per-label metrics
    label_metrics = evaluate_by_group(df_valid, 'label', ref_col,
                                      list(pred_cols.values()))

    # Per-phylum metrics (if column exists)
    if 'q_phylum' in df_valid.columns:
        phylum_metrics = evaluate_by_group(df_valid, 'q_phylum', ref_col,
                                           list(pred_cols.values()))
    else:
        phylum_metrics = {}

    report = {
        'total_pairs': len(df_valid),
        'global': global_metrics,
        'by_label': label_metrics,
        'by_phylum': phylum_metrics,
    }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(report, f, indent=2)
    print(f'\nSaved evaluation JSON: {out_path}')

    # Print summary
    print('\n=== Global metrics ===')
    for name, metrics in global_metrics.items():
        print(f"{name}: n={metrics['n']}, MAE={metrics['mae']*100:.3f}%, "
              f"RMSE={metrics['rmse']*100:.3f}%, r={metrics['pearson']:.4f}")

    print('\n=== Per-label MAE ===')
    for label in ['high', 'mid_high', 'mid', 'low']:
        if label in label_metrics:
            print(f"{label}:")
            for name, col in pred_cols.items():
                m = label_metrics[label].get(col, {})
                if m.get('mae') is not None:
                    print(f"  {name}: MAE={m['mae']*100:.3f}% (n={m['n']})")

    # Figures
    fig_dir = Path(args.figures)
    fig_dir.mkdir(parents=True, exist_ok=True)

    plot_scatter(df_valid, ref_col, pred_cols,
                 fig_dir / 'gtdb_r207_scatter.png')
    plot_error_by_label(df_valid, ref_col, pred_cols,
                        fig_dir / 'gtdb_r207_error_by_label.png')
    if 'q_phylum' in df_valid.columns:
        plot_phylum_error(df_valid, ref_col, pred_cols,
                          fig_dir / 'gtdb_r207_phylum_error.png',
                          min_pairs=args.min_phylum_pairs)


if __name__ == '__main__':
    main()
