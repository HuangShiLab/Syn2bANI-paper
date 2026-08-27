#!/usr/bin/env python3
"""Sample a stratified subset of pairs for FastANI reference computation.

The subset is stratified by skani ANI (or by the pair label if present) so that
the resulting FastANI reference covers the full ANI range and is suitable for
training the GBRT debiasing model.

Usage:
  python3 scripts/sample_fastani_subset.py \
    --matrix results/matrix.tsv \
    --output results/pairs_fastani_subset.tsv \
    --n-per-bin 200
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def sample_by_skani_ani(df: pd.DataFrame, n_per_bin: int, seed: int = 42) -> pd.DataFrame:
    """Stratified sample by skani ANI bins."""
    rng = np.random.default_rng(seed)
    # Define bins in fraction space (skani_ani is 0-1)
    bins = [0.0, 0.75, 0.80, 0.85, 0.90, 0.95, 0.975, 1.0]
    bin_labels = ['<0.75', '0.75-0.80', '0.80-0.85', '0.85-0.90',
                  '0.90-0.95', '0.95-0.975', '0.975-1.0']

    df = df[df['skani_ani'].notna()].copy()
    df['ani_bin'] = pd.cut(df['skani_ani'], bins=bins, labels=bin_labels,
                           include_lowest=True, right=False)

    sampled = []
    for label, group in df.groupby('ani_bin', observed=False):
        if len(group) == 0:
            continue
        n = min(n_per_bin, len(group))
        sampled.append(group.sample(n=n, random_state=rng.integers(0, 2**31)))

    if not sampled:
        return pd.DataFrame(columns=df.columns)
    return pd.concat(sampled, ignore_index=True)


def sample_by_label(df: pd.DataFrame, n_per_label: int, seed: int = 42) -> pd.DataFrame:
    """Stratified sample by existing relatedness label."""
    rng = np.random.default_rng(seed)
    if 'label' not in df.columns:
        return pd.DataFrame(columns=df.columns)

    sampled = []
    for label, group in df.groupby('label'):
        if len(group) == 0:
            continue
        n = min(n_per_label, len(group))
        sampled.append(group.sample(n=n, random_state=rng.integers(0, 2**31)))

    if not sampled:
        return pd.DataFrame(columns=df.columns)
    return pd.concat(sampled, ignore_index=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix', required=True,
                        help='Merged matrix with skani_ani')
    parser.add_argument('--output', required=True,
                        help='Output pair TSV for FastANI')
    parser.add_argument('--n-per-bin', type=int, default=200,
                        help='Pairs per skani ANI bin')
    parser.add_argument('--n-per-label', type=int, default=None,
                        help='Pairs per label (overrides --n-per-bin if set)')
    parser.add_argument('--seed', type=int, default=42)
    args = parser.parse_args()

    df = pd.read_csv(args.matrix, sep='\t', low_memory=False)
    print(f'Loaded {len(df)} rows from {args.matrix}')

    if args.n_per_label is not None and 'label' in df.columns:
        sampled = sample_by_label(df, args.n_per_label, args.seed)
    else:
        sampled = sample_by_skani_ani(df, args.n_per_bin, args.seed)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    sampled[['query', 'reference']].to_csv(out_path, sep='\t', index=False)
    print(f'Saved {len(sampled)} pairs to {out_path}')

    if 'ani_bin' in sampled.columns:
        print('\nSamples per bin:')
        print(sampled['ani_bin'].value_counts().sort_index())
    if 'label' in sampled.columns:
        print('\nSamples per label:')
        print(sampled['label'].value_counts())


if __name__ == '__main__':
    main()
