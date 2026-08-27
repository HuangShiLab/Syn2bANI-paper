#!/usr/bin/env python3
"""Generate GTDB-R207 benchmark figures."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RESULTS = Path('results')
FIGURES = Path('figures')
FIGURES.mkdir(exist_ok=True)


def load_eval():
    with open(RESULTS / 'evaluation_gtdb_r207_1k_optimized.json') as f:
        return json.load(f)


def load_matrix():
    return pd.read_csv(RESULTS / 'matrix_gtdb_r207_1k_optimized.tsv', sep='\t')


def plot_scatter(df):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    truth = df['fastani_ani'].values * 100
    methods = [
        ('Syn2bANI raw', df['s2b_raw_ani'].values * 100, '#E84855'),
        ('Syn2bANI GBRT', df['s2b_gbrt_ani'].values * 100, '#2E86AB'),
        ('skani', df['skani_ani'].values * 100, '#F6AE2D'),
    ]
    for ax, (name, pred, color) in zip(axes, methods):
        ax.scatter(truth, pred, c=color, alpha=0.4, s=25, edgecolors='none')
        lim = [min(truth.min(), pred.min()) - 1, max(truth.max(), pred.max()) + 1]
        ax.plot([0, 100], [0, 100], 'k--', lw=1)
        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xlabel('FastANI ANI (%)', fontsize=11)
        ax.set_ylabel(f'{name} ANI (%)', fontsize=11)
        ax.set_title(name, fontsize=12)
        mae = np.mean(np.abs(truth - pred))
        ax.text(0.05, 0.95, f'MAE = {mae:.2f}%', transform=ax.transAxes,
                fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
    plt.tight_layout()
    fig.savefig(FIGURES / 'gtdb_r207_scatter.png', dpi=300)
    plt.close(fig)
    print(f"Saved {FIGURES / 'gtdb_r207_scatter.png'}")


def plot_error_by_label(eval_data):
    labels = ['high', 'mid_high', 'mid', 'low']
    label_names = ['High (≥95%)', 'Mid-high (90–95%)', 'Mid (85–90%)', 'Low (<85%)']
    raw = [eval_data['by_label'][l]['s2b_raw_ani']['mae'] * 100 for l in labels]
    gbrt = [eval_data['by_label'][l]['s2b_gbrt_ani']['mae'] * 100 for l in labels]
    skani = [eval_data['by_label'][l]['skani_ani']['mae'] * 100 for l in labels]

    x = np.arange(len(labels))
    width = 0.25
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(x - width, raw, width, label='Syn2bANI raw', color='#E84855')
    ax.bar(x, gbrt, width, label='Syn2bANI GBRT', color='#2E86AB')
    ax.bar(x + width, skani, width, label='skani', color='#F6AE2D')
    ax.set_ylabel('MAE vs FastANI (%)', fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(label_names)
    ax.set_yscale('log')
    ax.set_ylim(0.001, 100)
    ax.legend()
    ax.set_title('GTDB-R207: ANI error by relatedness label', fontsize=13)
    plt.tight_layout()
    fig.savefig(FIGURES / 'gtdb_r207_error_by_label.png', dpi=300)
    plt.close(fig)
    print(f"Saved {FIGURES / 'gtdb_r207_error_by_label.png'}")


def plot_phylum_error(eval_data, top_n=15):
    by_phylum = eval_data['by_phylum']
    rows = []
    for phylum, vals in by_phylum.items():
        n = vals['s2b_raw_ani']['n']
        if n >= 3:
            rows.append({
                'phylum': phylum.replace('p__', ''),
                'n': n,
                'raw_mae': vals['s2b_raw_ani']['mae'] * 100,
                'gbrt_mae': vals['s2b_gbrt_ani']['mae'] * 100,
                'skani_mae': vals['skani_ani']['mae'] * 100,
            })
    df = pd.DataFrame(rows).sort_values('n', ascending=False).head(top_n)
    df = df.sort_values('gbrt_mae', ascending=True)

    fig, ax = plt.subplots(figsize=(10, 7))
    y = np.arange(len(df))
    ax.barh(y - 0.25, df['raw_mae'], 0.25, label='Syn2bANI raw', color='#E84855')
    ax.barh(y, df['gbrt_mae'], 0.25, label='Syn2bANI GBRT', color='#2E86AB')
    ax.barh(y + 0.25, df['skani_mae'], 0.25, label='skani', color='#F6AE2D')
    ax.set_yticks(y)
    ax.set_yticklabels([f"{p} (n={n})" for p, n in zip(df['phylum'], df['n'])])
    ax.set_xlabel('MAE vs FastANI (%)', fontsize=12)
    ax.set_title(f'GTDB-R207: per-phylum error (top {top_n} by sample size, n≥3)', fontsize=13)
    ax.legend()
    ax.set_xlim(0, max(df[['raw_mae', 'gbrt_mae', 'skani_mae']].max()) * 1.15)
    plt.tight_layout()
    fig.savefig(FIGURES / 'gtdb_r207_phylum_error.png', dpi=300)
    plt.close(fig)
    print(f"Saved {FIGURES / 'gtdb_r207_phylum_error.png'}")


def main():
    eval_data = load_eval()
    df = load_matrix()
    plot_scatter(df)
    plot_error_by_label(eval_data)
    plot_phylum_error(eval_data)


if __name__ == '__main__':
    main()
