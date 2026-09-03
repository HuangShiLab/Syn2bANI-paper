#!/usr/bin/env python3
"""Diagnose the closed-genome inverted-fraction distribution and seed-pair overlap."""

import argparse
import sys
from pathlib import Path

import pandas as pd


def load_closed_tsv(path):
    df = pd.read_csv(path, sep='\t')
    # pairid should be acc1__acc2
    df[['q_acc', 'r_acc']] = df['pairid'].str.split('__', n=1, expand=True)
    return df


def diagnose(seed_pairs_path, closed_outputs_path, genomes_path, out_prefix):
    seed = pd.read_csv(seed_pairs_path, sep='\t')
    out = load_closed_tsv(closed_outputs_path)
    genomes = pd.read_csv(genomes_path, sep='\t')

    seed_ids = set(seed['pairid'])
    out_ids = set(out['pairid'])

    print(f'Seed pairs: {len(seed_ids)}')
    print(f'Closed outputs: {len(out)}')
    print(f'Seed IDs present in output: {len(seed_ids & out_ids)}')
    print(f'Seed IDs missing from output: {len(seed_ids - out_ids)}')

    if seed_ids - out_ids:
        print('\nFirst 20 missing seed pair IDs:')
        for pid in sorted(seed_ids - out_ids)[:20]:
            print(f'  {pid}')

    # Check accession-level overlap
    seed_accs = set(seed['q_acc']) | set(seed['r_acc'])
    out_accs = set(out['q_acc']) | set(out['r_acc'])
    print(f'\nSeed accessions: {len(seed_accs)}')
    print(f'Output accessions: {len(out_accs)}')
    print(f'Seed accessions in output: {len(seed_accs & out_accs)}')
    print(f'Seed accessions missing from output: {len(seed_accs - out_accs)}')

    if seed_accs - out_accs:
        print('\nFirst 20 missing seed accessions:')
        for acc in sorted(seed_accs - out_accs)[:20]:
            print(f'  {acc}')

    # Annotate outputs with genome metadata
    acc_to_species = dict(zip(genomes['acc'], genomes['species']))
    acc_to_contigs = dict(zip(genomes['acc'], genomes['contig_count']))
    acc_to_level = dict(zip(genomes['acc'], genomes.get('ncbi_assembly_level', [])))

    out['q_species'] = out['q_acc'].map(acc_to_species)
    out['r_species'] = out['r_acc'].map(acc_to_species)
    out['q_contigs'] = out['q_acc'].map(acc_to_contigs)
    out['r_contigs'] = out['r_acc'].map(acc_to_contigs)
    out['same_species'] = out['q_species'] == out['r_species']

    # Distribution by species
    print('\n=== Inverted fraction by species (top 15 by median) ===')
    species_stats = out.groupby('q_species')['syn2b_raw_inverted_fraction'].agg(
        n='size', median='median', mean='mean', q90=lambda x: x.quantile(0.90)
    ).sort_values('median', ascending=False).head(15)
    print(species_stats.to_string())

    # P. aeruginosa detail
    pa = out[out['q_species'] == 'Pseudomonas aeruginosa']
    print(f"\nPseudomonas aeruginosa pairs: {len(pa)}")
    if len(pa):
        print(pa['syn2b_raw_inverted_fraction'].describe())

    # Correlation with contig count
    if out['q_contigs'].notna().any():
        corr = out[['syn2b_raw_inverted_fraction', 'q_contigs', 'r_contigs']].corr()
        print('\n=== Correlation with contig count ===')
        print(corr.to_string())

    # Compare raw vs corrected inverted fraction
    diff = (out['syn2b_raw_inverted_fraction'] - out['syn2b_inverted_fraction']).abs()
    print(f"\nraw vs corrected inverted_fraction differ in { (diff > 1e-6).sum() } / {len(out)} pairs")
    print(f"Max absolute difference: {diff.max():.4f}")

    # Save diagnostic tables
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    missing = seed[~seed['pairid'].isin(out_ids)].copy()
    missing.to_csv(out_prefix.with_suffix('.missing_seed_pairs.tsv'), sep='\t', index=False)

    species_stats.to_csv(out_prefix.with_suffix('.species_stats.tsv'), sep='\t')

    print(f"\nWrote diagnostic tables to {out_prefix}.*")


def main():
    parser = argparse.ArgumentParser(description='Diagnose closed-genome inversion outputs')
    parser.add_argument('--seed-pairs', required=True, type=Path)
    parser.add_argument('--closed-outputs', required=True, type=Path)
    parser.add_argument('--genomes', required=True, type=Path)
    parser.add_argument('--out-prefix', default='results/gtdb50k/closed_inversion_diagnostic', type=Path)
    args = parser.parse_args()
    diagnose(args.seed_pairs, args.closed_outputs, args.genomes, args.out_prefix)


if __name__ == '__main__':
    main()
