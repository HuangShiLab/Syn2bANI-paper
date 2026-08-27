#!/usr/bin/env python3
"""Analyze single-enzyme sweep and simulate multi-enzyme combinations."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def mash_from_af(af_q, af_r, tag_len):
    c1 = np.maximum(af_q, 1e-10)
    c2 = np.maximum(af_r, 1e-10)
    geo = np.sqrt(c1 * c2)
    return np.clip(1.0 + np.log(geo) / tag_len, 0.0, 1.0)


# Effective tag lengths from registry.rs
TAG_LENGTHS = {
    'BcgI': 32, 'AlfI': 32, 'AloI': 27, 'BaeI': 28, 'BplI': 27, 'BsaXI': 27,
    'BslFI': 21, 'Bsp24I': 27, 'CjeI': 28, 'CjePI': 27, 'CspCI': 33, 'FalI': 27,
    'HaeIV': 27, 'Hin4I': 27, 'PpiI': 28, 'PsrI': 27,
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sweep', required=True)
    parser.add_argument('--pairs', required=True)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.sweep, sep='\t', low_memory=False)
    if 'error' in df.columns:
        df = df[df['error'].isna()].copy()
    df['tag_len'] = df['enzyme'].map(TAG_LENGTHS)

    pairs = pd.read_csv(args.pairs, sep='\t', low_memory=False)
    pairs = pairs[['query', 'reference', 'fastani_ani', 'skani_ani']].copy()

    # Per-enzyme summary
    summary_rows = []
    for enzyme in sorted(df['enzyme'].unique()):
        sub = df[df['enzyme'] == enzyme]
        merged = sub.merge(pairs, on=['query', 'reference'])
        if len(merged) == 0:
            continue
        tag_len = TAG_LENGTHS.get(enzyme, 27)
        # Recompute mash_ani with correct tag length just in case
        recomp_mash = mash_from_af(merged['af_q'], merged['af_r'], tag_len)
        row = {
            'enzyme': enzyme,
            'n': len(merged),
            'mean_shared_tags': merged['shared_tags'].mean(),
            'mean_af_min': np.minimum(merged['af_q'], merged['af_r']).mean(),
            'mean_raw_err_pct': (merged['raw_ani'] - merged['fastani_ani']).mean() * 100,
            'mae_raw_pct': np.abs(merged['raw_ani'] - merged['fastani_ani']).mean() * 100,
            'mae_mash_pct': np.abs(recomp_mash - merged['fastani_ani']).mean() * 100,
            'mae_corrected_pct': np.abs(merged['corrected_ani'] - merged['fastani_ani']).mean() * 100,
        }
        summary_rows.append(row)
    summary = pd.DataFrame(summary_rows).sort_values('mae_mash_pct')
    summary_path = outdir / 'enzyme_summary.tsv'
    summary.to_csv(summary_path, sep='\t', index=False, float_format='%.4f')
    print(f'Saved: {summary_path}')
    print(summary.to_string(index=False, float_format='%.3f'))

    # Simulate multi-enzyme combinations by summing shared_tags
    # For each pair, sort enzymes by shared_tags contribution and compute cumulative metrics
    sim_records = []
    for (q, r), group in df.groupby(['query', 'reference']):
        if len(group) == 0:
            continue
        ref = pairs[(pairs['query'] == q) & (pairs['reference'] == r)]
        if len(ref) == 0:
            continue
        fastani = ref['fastani_ani'].values[0]
        skani = ref['skani_ani'].values[0]

        # Sort by shared_tags descending
        sorted_enz = group.sort_values('shared_tags', ascending=False).reset_index(drop=True)

        cum_shared = 0
        # For AF, we approximate total tags. Note: this is an upper-bound approximation
        # because unmatched sets from different enzymes may overlap.
        for k in range(1, min(5, len(sorted_enz)) + 1):
            sub = sorted_enz.iloc[:k]
            shared_sum = sub['shared_tags'].sum()
            # Approximate combined af as shared_sum / (shared_sum + sum of unmatched)
            unmatched_q = (sub['shared_tags'] / np.maximum(sub['af_q'], 1e-10) - sub['shared_tags']).sum()
            unmatched_r = (sub['shared_tags'] / np.maximum(sub['af_r'], 1e-10) - sub['shared_tags']).sum()
            af_q = shared_sum / max(shared_sum + unmatched_q, 1e-10)
            af_r = shared_sum / max(shared_sum + unmatched_r, 1e-10)
            # Use average tag length of selected enzymes
            avg_tag_len = sub['tag_len'].mean()
            sim_mash = mash_from_af(af_q, af_r, avg_tag_len)
            sim_records.append({
                'query': q,
                'reference': r,
                'fastani_ani': fastani,
                'skani_ani': skani,
                'n_enzymes': k,
                'enzymes': ','.join(sub['enzyme'].tolist()),
                'shared_tags': shared_sum,
                'af_min': min(af_q, af_r),
                'sim_mash_ani': sim_mash,
            })
    sim = pd.DataFrame(sim_records)
    sim_path = outdir / 'simulated_combinations.tsv'
    sim.to_csv(sim_path, sep='\t', index=False, float_format='%.6f')

    # Summary by n_enzymes
    combo_summary = sim.groupby('n_enzymes').agg(
        n=('query', 'count'),
        mean_shared=('shared_tags', 'mean'),
        mean_af_min=('af_min', 'mean'),
        mae_sim_mash=('sim_mash_ani', lambda x: np.abs(x - sim.loc[x.index, 'fastani_ani']).mean() * 100),
        mean_sim_mash_err=('sim_mash_ani', lambda x: (x - sim.loc[x.index, 'fastani_ani']).mean() * 100),
    ).reset_index()
    combo_summary_path = outdir / 'combination_summary.tsv'
    combo_summary.to_csv(combo_summary_path, sep='\t', index=False, float_format='%.4f')
    print(f'\nSaved: {combo_summary_path}')
    print(combo_summary.to_string(index=False, float_format='.3f'))

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].bar(summary['enzyme'], summary['mean_shared_tags'])
    axes[0].set_ylabel('Mean shared tags')
    axes[0].set_title('Shared tags per enzyme (validation pairs)')
    axes[0].tick_params(axis='x', rotation=45)

    axes[1].bar(summary['enzyme'], summary['mae_mash_pct'])
    axes[1].set_ylabel('MAE vs FastANI (%)')
    axes[1].set_title('Mash ANI MAE per enzyme')
    axes[1].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    fig.savefig(outdir / 'enzyme_sweep_summary.png', dpi=150)
    plt.close(fig)
    print(f'Saved plot: {outdir / "enzyme_sweep_summary.png"}')


if __name__ == '__main__':
    main()
