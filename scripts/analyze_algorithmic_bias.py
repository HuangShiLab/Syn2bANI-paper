#!/usr/bin/env python3
"""Analyze algorithmic sources of Syn2bANI mid-ANI bias vs FastANI/skani."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load_train(path):
    df = pd.read_csv(path, sep='\t', low_memory=False)
    df = df[df['fastani_ani'].notna()].copy()
    return df


def load_validation(path):
    df = pd.read_csv(path, sep='\t', low_memory=False)
    df = df[df['fastani_ani'].notna()].copy()
    return df


def add_derived(df, prefix):
    df = df.copy()
    raw = df[f'{prefix}_raw_ani']
    mash = df[f'{prefix}_mash_ani']
    af_q = df[f'{prefix}_af_q']
    af_r = df[f'{prefix}_af_r']
    shared = df[f'{prefix}_shared_tags']

    df[f'{prefix}_raw_err'] = (raw - df['fastani_ani']) * 100
    df[f'{prefix}_mash_err'] = (mash - df['fastani_ani']) * 100
    df[f'{prefix}_af_min'] = np.minimum(af_q, af_r)
    df[f'{prefix}_af_max'] = np.maximum(af_q, af_r)
    df[f'{prefix}_af_geo'] = np.sqrt(np.maximum(af_q, 1e-10) * np.maximum(af_r, 1e-10))
    # Estimated total tags and unmatched tags from shared tags / af
    df[f'{prefix}_total_tags_q'] = shared / np.maximum(af_q, 1e-10)
    df[f'{prefix}_total_tags_r'] = shared / np.maximum(af_r, 1e-10)
    df[f'{prefix}_unmatched_q'] = df[f'{prefix}_total_tags_q'] - shared
    df[f'{prefix}_unmatched_r'] = df[f'{prefix}_total_tags_r'] - shared
    df[f'{prefix}_unmatched_ratio'] = (df[f'{prefix}_unmatched_q'] + df[f'{prefix}_unmatched_r']) / (df[f'{prefix}_total_tags_q'] + df[f'{prefix}_total_tags_r'])
    return df


def bin_summary(df, prefix, ani_col='fastani_ani'):
    bins = [0.75, 0.80, 0.85, 0.88, 0.90, 0.92, 0.95, 1.00]
    labels = ['<80', '80-85', '85-88', '88-90', '90-92', '92-95', '95-100']
    df = df.copy()
    df['ani_bin'] = pd.cut(df[ani_col], bins=bins, labels=labels, right=False)
    summary = df.groupby('ani_bin').agg(
        n=(ani_col, 'count'),
        mean_fastani_ani=(ani_col, 'mean'),
        mean_raw_err=(f'{prefix}_raw_err', 'mean'),
        mae_raw=(f'{prefix}_raw_err', lambda x: x.abs().mean()),
        mean_mash_err=(f'{prefix}_mash_err', 'mean'),
        mae_mash=(f'{prefix}_mash_err', lambda x: x.abs().mean()),
        mean_shared=(f'{prefix}_shared_tags', 'mean'),
        median_af_min=(f'{prefix}_af_min', 'median'),
        mean_af_min=(f'{prefix}_af_min', 'mean'),
        mean_unmatched_ratio=(f'{prefix}_unmatched_ratio', 'mean'),
    ).reset_index()
    return summary


def analyze_matrix_100k(path):
    df = pd.read_csv(path, sep='\t', low_memory=False)
    df = df.dropna(subset=['fastani_ani']).copy()
    df['ani_pct'] = df['fastani_ani'] * 100
    df['s2b_af_min'] = np.minimum(df['s2b_af_q'], df['s2b_af_r'])

    bins = [0, 75, 80, 85, 90, 95, 100]
    df['bin'] = pd.cut(df['ani_pct'], bins=bins)
    summary = df.groupby('bin').agg(
        n=('fastani_ani', 'count'),
        mean_ani=('ani_pct', 'mean'),
        median_skani_af=('skani_align_frac', 'median'),
        mean_skani_af=('skani_align_frac', 'mean'),
        median_s2b_af_min=('s2b_af_min', 'median'),
        mean_s2b_af_min=('s2b_af_min', 'mean'),
        median_s2b_shared=('s2b_shared_tags', 'median'),
    ).reset_index()

    low_af = df[df['skani_align_frac'] < 0.5]
    low_af_s2b = df[df['s2b_af_min'] < 0.5]

    thresholds = []
    sorted_df = df.sort_values('ani_pct')
    for threshold in [95, 90, 85, 80, 75, 70]:
        sub = sorted_df[sorted_df['ani_pct'] <= threshold]
        if len(sub) > 0:
            pct_skani = (sub['skani_align_frac'] < 0.5).mean() * 100
            pct_s2b = (sub['s2b_af_min'] < 0.5).mean() * 100
            thresholds.append({
                'threshold': threshold,
                'pct_skani_af_below_50': pct_skani,
                'pct_s2b_af_min_below_50': pct_s2b,
            })

    return {
        'summary': summary,
        'low_af': low_af,
        'low_af_s2b': low_af_s2b,
        'thresholds': pd.DataFrame(thresholds),
        'df': df,
    }


def per_genus_summary(df, prefix, min_n=2):
    rows = []
    for genus in sorted(df['q_genus'].unique()):
        sub = df[df['q_genus'] == genus]
        if len(sub) < min_n:
            continue
        row = {
            'genus': genus,
            'n': len(sub),
            'mean_fastani_ani': sub['fastani_ani'].mean(),
            'mean_raw_err': sub[f'{prefix}_raw_err'].mean(),
            'mae_raw': sub[f'{prefix}_raw_err'].abs().mean(),
            'mean_mash_err': sub[f'{prefix}_mash_err'].mean(),
            'mae_mash': sub[f'{prefix}_mash_err'].abs().mean(),
            'mean_shared': sub[f'{prefix}_shared_tags'].mean(),
            'mean_af_min': sub[f'{prefix}_af_min'].mean(),
        }
        rows.append(row)
    return pd.DataFrame(rows).sort_values('mae_raw', ascending=False)


def plot_scatter(df, prefix, outdir):
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    axes = axes.flatten()

    x = df['fastani_ani'] * 100
    axes[0].scatter(x, df[f'{prefix}_raw_err'], alpha=0.4, s=8)
    axes[0].axhline(0, color='red', linestyle='--')
    axes[0].set_xlabel('FastANI ANI (%)')
    axes[0].set_ylabel('Raw ANI error (%)')
    axes[0].set_title(f'{prefix}: raw ANI error vs FastANI')

    axes[1].scatter(x, df[f'{prefix}_mash_err'], alpha=0.4, s=8)
    axes[1].axhline(0, color='red', linestyle='--')
    axes[1].set_xlabel('FastANI ANI (%)')
    axes[1].set_ylabel('Mash ANI error (%)')
    axes[1].set_title(f'{prefix}: mash ANI error vs FastANI')

    axes[2].scatter(x, df[f'{prefix}_shared_tags'], alpha=0.4, s=8)
    axes[2].set_xlabel('FastANI ANI (%)')
    axes[2].set_ylabel('Shared tags')
    axes[2].set_title(f'{prefix}: shared tags vs FastANI')

    axes[3].scatter(x, df[f'{prefix}_af_min'] * 100, alpha=0.4, s=8)
    axes[3].axhline(50, color='red', linestyle='--')
    axes[3].set_xlabel('FastANI ANI (%)')
    axes[3].set_ylabel('Min alignment fraction (%)')
    axes[3].set_title(f'{prefix}: min AF vs FastANI')

    axes[4].scatter(df[f'{prefix}_af_min'] * 100, df[f'{prefix}_raw_err'], alpha=0.4, s=8)
    axes[4].axhline(0, color='red', linestyle='--')
    axes[4].set_xlabel('Min alignment fraction (%)')
    axes[4].set_ylabel('Raw ANI error (%)')
    axes[4].set_title(f'{prefix}: raw error vs min AF')

    axes[5].scatter(df[f'{prefix}_unmatched_ratio'] * 100, df[f'{prefix}_raw_err'], alpha=0.4, s=8)
    axes[5].axhline(0, color='red', linestyle='--')
    axes[5].set_xlabel('Unmatched tags ratio (%)')
    axes[5].set_ylabel('Raw ANI error (%)')
    axes[5].set_title(f'{prefix}: raw error vs unmatched ratio')

    plt.tight_layout()
    outpath = Path(outdir) / f'{prefix}_bias_diagnostics.png'
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    return outpath


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--train', default='/lustre1/g/aos_shihuang/data/gtdb-r207/train_combined.tsv')
    parser.add_argument('--val', default='/lustre1/g/aos_shihuang/data/validation_mid_ani/val_combined.tsv')
    parser.add_argument('--matrix100k', default='/lustre1/g/aos_shihuang/Syn2bANI-paper/results/matrix_gtdb_r207_100k.tsv')
    parser.add_argument('--outdir', default='/lustre1/g/aos_shihuang/Syn2bANI-paper/results/algorithm_analysis')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    train = load_train(args.train)
    val = load_validation(args.val)

    report_lines = []
    report_lines.append('# Syn2bANI Algorithmic Bias Analysis')
    report_lines.append('')

    # 1. 100k matrix alignment fraction analysis
    report_lines.append('## 1. Alignment fraction vs ANI (GTDB-R207 100k matrix)')
    report_lines.append('')
    matrix = analyze_matrix_100k(args.matrix100k)
    report_lines.append('### ANI bins')
    report_lines.append(matrix['summary'].to_string(index=False, float_format='%.3f'))
    report_lines.append('')
    report_lines.append(f"Pairs with skani align frac < 50%: {len(matrix['low_af'])} / {len(matrix['df'])} ({len(matrix['low_af'])/len(matrix['df'])*100:.1f}%)")
    if len(matrix['low_af']) > 0:
        report_lines.append(f"  ANI range: {matrix['low_af']['ani_pct'].min():.2f}% - {matrix['low_af']['ani_pct'].max():.2f}%, mean: {matrix['low_af']['ani_pct'].mean():.2f}%")
    report_lines.append(f"Pairs with Syn2bANI min AF < 50%: {len(matrix['low_af_s2b'])} / {len(matrix['df'])} ({len(matrix['low_af_s2b'])/len(matrix['df'])*100:.1f}%)")
    if len(matrix['low_af_s2b']) > 0:
        report_lines.append(f"  ANI range: {matrix['low_af_s2b']['ani_pct'].min():.2f}% - {matrix['low_af_s2b']['ani_pct'].max():.2f}%, mean: {matrix['low_af_s2b']['ani_pct'].mean():.2f}%")
    report_lines.append('')
    report_lines.append('### Thresholds: fraction of pairs with AF < 50% below given ANI')
    report_lines.append(matrix['thresholds'].to_string(index=False, float_format='%.1f'))
    report_lines.append('')

    # 2. BcgI / CjePI comparison
    for prefix in ['bcgi', 'cjepi']:
        report_lines.append(f'## 2. {prefix.upper()} training analysis')
        report_lines.append('')
        train_d = add_derived(train, prefix)
        val_d = add_derived(val, prefix)

        report_lines.append('### Overall error')
        for label, df_d in [('train', train_d), ('validation', val_d)]:
            raw_mae = df_d[f'{prefix}_raw_err'].abs().mean()
            mash_mae = df_d[f'{prefix}_mash_err'].abs().mean()
            raw_mean_err = df_d[f'{prefix}_raw_err'].mean()
            mash_mean_err = df_d[f'{prefix}_mash_err'].mean()
            report_lines.append(f"- {label}: raw MAE={raw_mae:.3f}%, mean raw err={raw_mean_err:+.3f}%; mash MAE={mash_mae:.3f}%, mean mash err={mash_mean_err:+.3f}%")
        report_lines.append('')

        report_lines.append('### By ANI bin (validation)')
        val_summary = bin_summary(val_d, prefix)
        report_lines.append(val_summary.to_string(index=False, float_format='%.3f'))
        report_lines.append('')

        report_lines.append('### By ANI bin (training)')
        train_summary = bin_summary(train_d, prefix)
        report_lines.append(train_summary.to_string(index=False, float_format='%.3f'))
        report_lines.append('')

        report_lines.append('### Per-genus validation summary')
        genus_summary = per_genus_summary(val_d, prefix, min_n=1)
        report_lines.append(genus_summary.to_string(index=False, float_format='%.3f'))
        report_lines.append('')

        # Save plots
        if len(train_d) > 0:
            plot_path = plot_scatter(train_d, prefix, outdir)
            report_lines.append(f"Scatter plots saved: {plot_path}")
            report_lines.append('')

        # Correlation between raw error and various quantities
        report_lines.append('### Correlations with raw error (validation)')
        corr_cols = [f'{prefix}_shared_tags', f'{prefix}_af_min', f'{prefix}_unmatched_ratio', f'{prefix}_mash_ani']
        for col in corr_cols:
            if col in val_d.columns:
                corr = val_d[f'{prefix}_raw_err'].corr(val_d[col])
                report_lines.append(f"- {col}: {corr:.3f}")
        report_lines.append('')

    # 3. Mash-like estimator formula check
    report_lines.append('## 3. Mash-like estimator formula check')
    report_lines.append('')
    report_lines.append('Mash ANI computed from geometric mean of AF: `1 + log(sqrt(af_q * af_r)) / tag_len`')
    report_lines.append('')
    for prefix in ['bcgi', 'cjepi']:
        train_d = add_derived(train, prefix)
        af_geo = np.sqrt(np.maximum(train_d[f'{prefix}_af_q'], 1e-10) * np.maximum(train_d[f'{prefix}_af_r'], 1e-10))
        tag_len = 32.0
        recomputed = 1.0 + np.log(af_geo) / tag_len
        max_diff = (recomputed - train_d[f'{prefix}_mash_ani']).abs().max()
        report_lines.append(f"- {prefix}: max recompute diff vs stored mash_ani = {max_diff:.6f}")
    report_lines.append('')

    # 4. Unmatched tags hypothesis
    report_lines.append('## 4. Unmatched tags hypothesis')
    report_lines.append('')
    report_lines.append('Estimated total tags = shared_tags / af; unmatched = total - shared.')
    report_lines.append('')
    for prefix in ['bcgi', 'cjepi']:
        train_d = add_derived(train, prefix)
        report_lines.append(f"### {prefix.upper()} (training)")
        report_lines.append(f"- mean shared tags: {train_d[f'{prefix}_shared_tags'].mean():.1f}")
        report_lines.append(f"- mean estimated total tags (query): {train_d[f'{prefix}_total_tags_q'].mean():.1f}")
        report_lines.append(f"- mean unmatched ratio: {train_d[f'{prefix}_unmatched_ratio'].mean():.3f}")
        report_lines.append(f"- correlation raw_err vs unmatched_ratio: {train_d[f'{prefix}_raw_err'].corr(train_d[f'{prefix}_unmatched_ratio']):.3f}")
        report_lines.append('')

    report_text = '\n'.join(report_lines)
    report_path = outdir / 'algorithm_bias_report.md'
    report_path.write_text(report_text)
    print(f'Report saved: {report_path}')

    # Print to stdout as well
    print(report_text)


if __name__ == '__main__':
    main()
