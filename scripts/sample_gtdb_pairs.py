#!/usr/bin/env python3
"""
sample_gtdb_pairs.py
Stratified sampling of genome pairs from GTDB-R207 for benchmarking.

Produces pairs.tsv with columns:
  query, reference, level (intra_species|intra_genus|intra_family|random)

Usage:
  python3 sample_gtdb_pairs.py \
    --bac-metadata metadata/bac120_metadata_r207.tsv \
    --ar-metadata metadata/ar53_metadata_r207.tsv \
    --output pairs_20k.tsv \
    --n-per-level 5000 \
    --seed 42
"""

import argparse
import random
import pandas as pd
from collections import defaultdict
from pathlib import Path


def load_metadata(bac_path: str, ar_path: str) -> pd.DataFrame:
    """Load and merge bacterial + archaeal metadata."""
    print(f"Loading bacterial metadata: {bac_path}")
    bac = pd.read_csv(bac_path, sep='\t', low_memory=False)
    print(f"  Rows: {len(bac)}")

    print(f"Loading archaeal metadata: {ar_path}")
    ar = pd.read_csv(ar_path, sep='\t', low_memory=False)
    print(f"  Rows: {len(ar)}")

    df = pd.concat([bac, ar], ignore_index=True)
    print(f"Combined: {len(df)} genomes")

    # Normalize accession column
    for col in ['ncbi_genbank_assembly_accession', 'accession',
                'genome_accession', 'gtdb_accession']:
        if col in df.columns:
            df['accession'] = df[col].astype(str)
            break

    # Keep relevant columns
    keep_cols = ['accession', 'gtdb_species', 'gtdb_genus', 'gtdb_family',
                 'gtdb_order', 'gtdb_class', 'gtdb_phylum',
                 'checkm_completeness', 'checkm_contamination',
                 'gc_percentage', 'genome_size']
    available = [c for c in keep_cols if c in df.columns]
    df = df[available].copy()

    # Quality filter: completeness > 90%, contamination < 5%
    if 'checkm_completeness' in df.columns:
        before = len(df)
        df = df[df['checkm_completeness'] > 90]
        print(f"Filtered by completeness >90%: {before} -> {len(df)}")

    if 'checkm_contamination' in df.columns:
        before = len(df)
        df = df[df['checkm_contamination'] < 5]
        print(f"Filtered by contamination <5%: {before} -> {len(df)}")

    df = df.dropna(subset=['accession', 'gtdb_species'])
    print(f"Final after cleaning: {len(df)} genomes")
    return df


def sample_intra_species(df: pd.DataFrame, n_target: int, rng: random.Random) -> list:
    """Sample pairs within the same species."""
    pairs = []
    groups = df.groupby('gtdb_species')
    species_list = [sp for sp, g in groups if len(g) >= 2]
    rng.shuffle(species_list)

    print(f"Species with >=2 genomes: {len(species_list)}")

    for sp in species_list:
        group = groups.get_group(sp)
        accs = group['accession'].tolist()
        rng.shuffle(accs)

        # Generate all pairwise combinations (capped)
        for i in range(len(accs)):
            for j in range(i + 1, len(accs)):
                pairs.append((accs[i], accs[j], 'intra_species'))
                if len([p for p in pairs if p[2] == 'intra_species']) >= n_target:
                    return pairs[:n_target]

    print(f"WARNING: Only found {len(pairs)} intra-species pairs (target: {n_target})")
    return pairs


def sample_intra_genus(df: pd.DataFrame, n_target: int, rng: random.Random) -> list:
    """Sample pairs within same genus, different species."""
    pairs = []
    genus_groups = df.groupby('gtdb_genus')
    genera = [g for g, grp in genus_groups
              if grp['gtdb_species'].nunique() >= 2]
    rng.shuffle(genera)

    print(f"Genera with >=2 species: {len(genera)}")

    for genus in genera:
        grp = genus_groups.get_group(genus)
        # Pick one representative per species
        reps = grp.groupby('gtdb_species')['accession'].apply(
            lambda x: rng.choice(x.tolist())
        ).tolist()
        rng.shuffle(reps)

        for i in range(min(len(reps) - 1, 100)):
            pairs.append((reps[i], reps[i + 1], 'intra_genus'))
            if len([p for p in pairs if p[2] == 'intra_genus']) >= n_target:
                return pairs[:n_target]

    print(f"WARNING: Only found {len(pairs)} intra-genus pairs (target: {n_target})")
    return pairs


def sample_intra_family(df: pd.DataFrame, n_target: int, rng: random.Random) -> list:
    """Sample pairs within same family, different genus."""
    pairs = []
    fam_groups = df.groupby('gtdb_family')
    families = [f for f, grp in fam_groups
                if grp['gtdb_genus'].nunique() >= 2]
    rng.shuffle(families)

    print(f"Families with >=2 genera: {len(families)}")

    for fam in families:
        grp = fam_groups.get_group(fam)
        reps = grp.groupby('gtdb_genus')['accession'].apply(
            lambda x: rng.choice(x.tolist())
        ).tolist()
        rng.shuffle(reps)

        for i in range(min(len(reps) - 1, 50)):
            pairs.append((reps[i], reps[i + 1], 'intra_family'))
            if len([p for p in pairs if p[2] == 'intra_family']) >= n_target:
                return pairs[:n_target]

    print(f"WARNING: Only found {len(pairs)} intra-family pairs (target: {n_target})")
    return pairs


def sample_random(df: pd.DataFrame, n_target: int, rng: random.Random) -> list:
    """Sample random cross-taxonomy pairs."""
    accs = df['accession'].tolist()
    pairs = []
    for _ in range(n_target * 3):  # oversample to avoid duplicates
        a, b = rng.sample(accs, 2)
        if a != b:
            pairs.append((a, b, 'random'))
        if len(pairs) >= n_target:
            break
    return pairs[:n_target]


def main():
    parser = argparse.ArgumentParser(
        description='Stratified sampling of GTDB-R207 genome pairs'
    )
    parser.add_argument('--bac-metadata', required=True,
                        help='Path to bac120_metadata_r207.tsv')
    parser.add_argument('--ar-metadata', required=True,
                        help='Path to ar53_metadata_r207.tsv')
    parser.add_argument('--output', required=True,
                        help='Output TSV file')
    parser.add_argument('--n-per-level', type=int, default=5000,
                        help='Target pairs per stratum (default: 5000)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed')
    args = parser.parse_args()

    rng = random.Random(args.seed)

    df = load_metadata(args.bac_metadata, args.ar_metadata)

    print(f"\nSampling {args.n_per_level} pairs per level (seed={args.seed})...")

    pairs = []
    pairs.extend(sample_intra_species(df, args.n_per_level, rng))
    pairs.extend(sample_intra_genus(df, args.n_per_level, rng))
    pairs.extend(sample_intra_family(df, args.n_per_level, rng))
    pairs.extend(sample_random(df, args.n_per_level, rng))

    out_df = pd.DataFrame(pairs, columns=['query', 'reference', 'level'])

    # Merge taxonomy info
    tax_cols = ['accession', 'gtdb_phylum', 'gtdb_class', 'gtdb_order',
                'gtdb_family', 'gtdb_genus', 'gtdb_species']
    tax_df = df[[c for c in tax_cols if c in df.columns]].copy()
    tax_df = tax_df.rename(columns={c: f"{c}_q" for c in tax_cols if c != 'accession'})

    out_df = out_df.merge(tax_df, left_on='query', right_on='accession', how='left')
    out_df = out_df.drop(columns=['accession'])

    tax_df = tax_df.rename(columns={c: c.replace('_q', '_r') for c in tax_df.columns if '_q' in c})
    tax_df = tax_df.rename(columns={'accession': 'reference'})
    out_df = out_df.merge(tax_df, on='reference', how='left')

    out_df.to_csv(args.output, sep='\t', index=False)

    print(f"\n{'='*50}")
    print(f"Output: {args.output}")
    print(f"Total pairs: {len(out_df)}")
    print(out_df['level'].value_counts())
    print(f"{'='*50}")
    print("\nNext: python3 scripts/run_benchmark_matrix.py --pairs", args.output)


if __name__ == '__main__':
    main()
