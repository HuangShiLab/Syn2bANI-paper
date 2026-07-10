#!/usr/bin/env python3
"""
evaluate_by_taxonomy.py
Layered error analysis by taxonomy + generate publication-ready figures.

Usage:
  python3 evaluate_by_taxonomy.py \
    --matrix results/matrix.tsv \
    --output figures/ \
    --report results/layered_error_report.txt
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import mean_absolute_error, mean_squared_error


# Set publication style
sns.set_style('ticks')
sns.set_context('paper', font_scale=1.3)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.family'] = 'sans-serif'


def load_matrix(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t', low_memory=False)
    # Compute errors
    df['error_raw'] = df['skani_ani'] - df['s2b_raw_ani']
    if 's2b_gbrt_ani' in df.columns and df['s2b_gbrt_ani'].notna().any():
        df['error_gbrt'] = df['skani_ani'] - df['s2b_gbrt_ani']
    else:
        df['error_gbrt'] = np.nan
    return df


def plot_scatter(df: pd.DataFrame, outdir: Path):
    """Figure 1: Scatter plot skani vs Syn2bANI (raw + GBRT)."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Raw
    ax = axes[0]
    sns.scatterplot(data=df, x='skani_ani', y='s2b_raw_ani',
                    hue='level', alpha=0.5, ax=ax, s=15)
    ax.plot([75, 100], [75, 100], 'k--', lw=1)
    ax.set_xlabel('skani ANI (%)')
    ax.set_ylabel('Syn2bANI raw (%)')
    ax.set_title('Raw Syn2bANI')
    ax.set_xlim(75, 100)
    ax.set_ylim(75, 100)

    # GBRT
    ax = axes[1]
    if df['error_gbrt'].notna().any():
        sns.scatterplot(data=df, x='skani_ani', y='s2b_gbrt_ani',
                        hue='level', alpha=0.5, ax=ax, s=15)
        ax.plot([75, 100], [75, 100], 'k--', lw=1)
        ax.set_xlabel('skani ANI (%)')
        ax.set_ylabel('Syn2bANI + GBRT (%)')
        ax.set_title('GBRT-corrected Syn2bANI')
        ax.set_xlim(75, 100)
        ax.set_ylim(75, 100)
    else:
        ax.text(0.5, 0.5, 'GBRT not available',
                ha='center', va='center', transform=ax.transAxes)

    plt.tight_layout()
    fig.savefig(outdir / 'fig1_scatter.png', bbox_inches='tight')
    fig.savefig(outdir / 'fig1_scatter.pdf', bbox_inches='tight')
    plt.close(fig)


def plot_error_distribution(df: pd.DataFrame, outdir: Path):
    """Figure 2: Error distribution by level."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Raw error
    ax = axes[0]
    sns.histplot(data=df, x='error_raw', hue='level',
                 bins=50, kde=True, ax=ax, alpha=0.6)
    ax.axvline(0, color='red', linestyle='--', lw=1)
    ax.set_xlabel('Error (skani - Syn2bANI raw, %)')
    ax.set_title('Raw Error Distribution')
    ax.set_xlim(-5, 5)

    # GBRT error
    ax = axes[1]
    if df['error_gbrt'].notna().any():
        sns.histplot(data=df, x='error_gbrt', hue='level',
                     bins=50, kde=True, ax=ax, alpha=0.6)
        ax.axvline(0, color='red', linestyle='--', lw=1)
        ax.set_xlabel('Error (skani - Syn2bANI GBRT, %)')
        ax.set_title('GBRT-corrected Error Distribution')
        ax.set_xlim(-5, 5)
    else:
        ax.text(0.5, 0.5, 'GBRT not available',
                ha='center', va='center', transform=ax.transAxes)

    plt.tight_layout()
    fig.savefig(outdir / 'fig2_error_dist.png', bbox_inches='tight')
    fig.savefig(outdir / 'fig2_error_dist.pdf', bbox_inches='tight')
    plt.close(fig)


def plot_by_taxonomy(df: pd.DataFrame, outdir: Path):
    """Figure 3: Boxplot of errors by phylum and class."""
    phylum_col_q = 'gtdb_phylum_q'
    if phylum_col_q not in df.columns:
        print(f"WARNING: {phylum_col_q} not in data. Skipping taxonomy plot.")
        return

    # Filter to phyla with enough samples
    phylum_counts = df[phylum_col_q].value_counts()
    top_phyla = phylum_counts[phylum_counts >= 50].index.tolist()
    df_sub = df[df[phylum_col_q].isin(top_phyla)].copy()

    fig, ax = plt.subplots(figsize=(14, 6))
    sns.boxplot(data=df_sub, x=phylum_col_q, y='error_raw', ax=ax)
    ax.axhline(0, color='red', linestyle='--', lw=1)
    ax.set_xlabel('Phylum')
    ax.set_ylabel('Error (skani - Syn2bANI raw, %)')
    ax.set_title('Error by Phylum (Raw Syn2bANI)')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    fig.savefig(outdir / 'fig3_error_by_phylum.png', bbox_inches='tight')
    fig.savefig(outdir / 'fig3_error_by_phylum.pdf', bbox_inches='tight')
    plt.close(fig)

    # GBRT version
    if df['error_gbrt'].notna().any():
        fig, ax = plt.subplots(figsize=(14, 6))
        sns.boxplot(data=df_sub, x=phylum_col_q, y='error_gbrt', ax=ax)
        ax.axhline(0, color='red', linestyle='--', lw=1)
        ax.set_xlabel('Phylum')
        ax.set_ylabel('Error (skani - Syn2bANI GBRT, %)')
        ax.set_title('Error by Phylum (GBRT-corrected)')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        fig.savefig(outdir / 'fig3_error_by_phylum_gbrt.png',
                    bbox_inches='tight')
        fig.savefig(outdir / 'fig3_error_by_phylum_gbrt.pdf',
                    bbox_inches='tight')
        plt.close(fig)


def plot_by_ani_range(df: pd.DataFrame, outdir: Path):
    """Figure 4: Error by ANI range."""
    df['ani_bin'] = pd.cut(df['skani_ani'],
                           bins=[0, 85, 90, 95, 97, 100],
                           labels=['<85', '85-90', '90-95', '95-97', '>97'])

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    sns.boxplot(data=df, x='ani_bin', y='error_raw', ax=axes[0])
    axes[0].axhline(0, color='red', linestyle='--', lw=1)
    axes[0].set_xlabel('skani ANI range (%)')
    axes[0].set_ylabel('Error (%)')
    axes[0].set_title('Raw Syn2bANI')

    if df['error_gbrt'].notna().any():
        sns.boxplot(data=df, x='ani_bin', y='error_gbrt', ax=axes[1])
        axes[1].axhline(0, color='red', linestyle='--', lw=1)
        axes[1].set_xlabel('skani ANI range (%)')
        axes[1].set_ylabel('Error (%)')
        axes[1].set_title('GBRT-corrected')

    plt.tight_layout()
    fig.savefig(outdir / 'fig4_error_by_ani_range.png', bbox_inches='tight')
    fig.savefig(outdir / 'fig4_error_by_ani_range.pdf', bbox_inches='tight')
    plt.close(fig)


def generate_report(df: pd.DataFrame, out_path: str):
    """Generate text report of layered errors."""
    lines = []
    lines.append("=" * 70)
    lines.append("Layered Error Analysis Report")
    lines.append("=" * 70)
    lines.append("")

    # Overall
    lines.append("--- Overall ---")
    mae_raw = mean_absolute_error(df['skani_ani'], df['s2b_raw_ani'])
    rmse_raw = np.sqrt(mean_squared_error(df['skani_ani'], df['s2b_raw_ani']))
    lines.append(f"Raw Syn2bANI MAE:  {mae_raw:.4f}%")
    lines.append(f"Raw Syn2bANI RMSE: {rmse_raw:.4f}%")

    if df['error_gbrt'].notna().any():
        mae_gbrt = mean_absolute_error(df['skani_ani'], df['s2b_gbrt_ani'])
        rmse_gbrt = np.sqrt(mean_squared_error(df['skani_ani'], df['s2b_gbrt_ani']))
        lines.append(f"GBRT MAE:          {mae_gbrt:.4f}%")
        lines.append(f"GBRT RMSE:         {rmse_gbrt:.4f}%")
        lines.append(f"MAE improvement:   {mae_raw - mae_gbrt:.4f}%")
    lines.append("")

    # By level
    lines.append("--- By Taxonomic Level ---")
    for level in df['level'].unique():
        subset = df[df['level'] == level]
        mae = mean_absolute_error(subset['skani_ani'], subset['s2b_raw_ani'])
        lines.append(f"  {level:20s}: n={len(subset):5d}, raw_MAE={mae:.4f}%")
        if df['error_gbrt'].notna().any():
            mae_g = mean_absolute_error(subset['skani_ani'], subset['s2b_gbrt_ani'])
            lines.append(f"  {'':20s}  gbrt_MAE={mae_g:.4f}%")
    lines.append("")

    # By phylum
    phylum_col = 'gtdb_phylum_q'
    if phylum_col in df.columns:
        lines.append("--- By Phylum (top 10) ---")
        phyla = df[phylum_col].value_counts().head(10).index
        for ph in phyla:
            subset = df[df[phylum_col] == ph]
            mae = mean_absolute_error(subset['skani_ani'], subset['s2b_raw_ani'])
            lines.append(f"  {ph:30s}: n={len(subset):4d}, raw_MAE={mae:.4f}%")
    lines.append("")

    # By ANI range
    lines.append("--- By ANI Range ---")
    df['ani_bin'] = pd.cut(df['skani_ani'],
                           bins=[0, 85, 90, 95, 97, 100],
                           labels=['<85', '85-90', '90-95', '95-97', '>97'])
    for b in df['ani_bin'].cat.categories:
        subset = df[df['ani_bin'] == b]
        if len(subset) == 0:
            continue
        mae = mean_absolute_error(subset['skani_ani'], subset['s2b_raw_ani'])
        lines.append(f"  {str(b):10s}: n={len(subset):5d}, raw_MAE={mae:.4f}%")

    lines.append("")
    lines.append("=" * 70)

    report = "\n".join(lines)
    with open(out_path, 'w') as f:
        f.write(report)
    print(report)
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix', required=True)
    parser.add_argument('--output', default='figures/')
    parser.add_argument('--report', default='results/layered_error_report.txt')
    args = parser.parse_args()

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading matrix: {args.matrix}")
    df = load_matrix(args.matrix)
    print(f"Samples: {len(df)}")

    print("\nGenerating figures...")
    plot_scatter(df, outdir)
    plot_error_distribution(df, outdir)
    plot_by_taxonomy(df, outdir)
    plot_by_ani_range(df, outdir)
    print(f"Figures saved to: {outdir}")

    print("\nGenerating report...")
    generate_report(df, args.report)
    print(f"Report saved to: {args.report}")


if __name__ == '__main__':
    main()
