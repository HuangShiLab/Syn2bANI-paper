#!/usr/bin/env python3
"""Generate validation scatter plots."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error, mean_squared_error


def plot_scatter(df, pred_col, ref_col, out_path, title):
    mask = df[ref_col].notna() & df[pred_col].notna()
    sub = df[mask]
    if len(sub) == 0:
        return

    mae = mean_absolute_error(sub[ref_col], sub[pred_col])
    rmse = np.sqrt(mean_squared_error(sub[ref_col], sub[pred_col]))
    r, _ = pearsonr(sub[ref_col], sub[pred_col])

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(sub[ref_col] * 100, sub[pred_col] * 100, alpha=0.6, edgecolors='none')
    lim = [min(sub[ref_col].min(), sub[pred_col].min()) * 100 - 1,
           max(sub[ref_col].max(), sub[pred_col].max()) * 100 + 1]
    ax.plot(lim, lim, 'k--', lw=1)
    ax.set_xlabel('FastANI ANI (%)')
    ax.set_ylabel(f'{pred_col} ANI (%)')
    ax.set_title(f'{title}\nMAE={mae*100:.2f}%, RMSE={rmse*100:.2f}%, r={r:.3f}')
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f'Saved: {out_path}')


def plot_residual(df, pred_col, ref_col, out_path, title):
    mask = df[ref_col].notna() & df[pred_col].notna()
    sub = df[mask]
    if len(sub) == 0:
        return

    residual = (sub[pred_col] - sub[ref_col]) * 100

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.scatter(sub[ref_col] * 100, residual, alpha=0.6, edgecolors='none')
    ax.axhline(0, color='k', linestyle='--', lw=1)
    ax.set_xlabel('FastANI ANI (%)')
    ax.set_ylabel(f'{pred_col} - FastANI (%)')
    ax.set_title(f'{title} Residuals')
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f'Saved: {out_path}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix', required=True)
    parser.add_argument('--out-dir', required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.matrix, sep='\t', low_memory=False)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    methods = {
        's2b_raw_ani': 'Syn2bANI raw',
        'gbrt_v4_ani': 'Syn2bANI GBRT v4',
        'skani_ani': 'skani',
        's2b_gbrt_ani': 'Syn2bANI built-in GBRT',
    }

    for col, name in methods.items():
        if col in df.columns:
            plot_scatter(df, col, 'fastani_ani',
                         out_dir / f'scatter_{col}.png', name)
            plot_residual(df, col, 'fastani_ani',
                          out_dir / f'residual_{col}.png', name)


if __name__ == '__main__':
    main()
