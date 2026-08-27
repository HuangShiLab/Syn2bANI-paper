#!/usr/bin/env python3
"""Run fastANI binary on all pairs in a matrix and merge results back."""
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Optional

import pandas as pd


def resolve_genome_path(accession: str, genomes_dir: Path) -> Optional[Path]:
    genomes_dir = Path(genomes_dir)
    direct = genomes_dir / f"{accession}.fna"
    if direct.exists():
        return direct
    matches = list(genomes_dir.glob(f"*{accession}*.fna"))
    if matches:
        return matches[0]
    return None


def run_fastani_pair(args: dict) -> dict:
    q, r = args['query'], args['reference']
    genomes_dir = args['genomes_dir']
    fastani_path = args['fastani_path']
    q_path = resolve_genome_path(q, genomes_dir)
    r_path = resolve_genome_path(r, genomes_dir)

    if q_path is None or r_path is None:
        return {'query': q, 'reference': r, 'fastani_ani': None,
                'error': 'missing_genome'}

    try:
        result = subprocess.run(
            [fastani_path, '-q', str(q_path), '-r', str(r_path), '-o', '/dev/stdout'],
            capture_output=True, text=True, timeout=300
        )
        lines = result.stdout.strip().split('\n')
        for line in lines:
            if not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) >= 3:
                try:
                    return {
                        'query': q, 'reference': r,
                        'fastani_ani': float(parts[2]) / 100.0,
                        'error': None
                    }
                except ValueError:
                    continue
        return {'query': q, 'reference': r, 'fastani_ani': None,
                'error': 'no_alignment'}
    except Exception as e:
        return {'query': q, 'reference': r, 'fastani_ani': None,
                'error': str(e)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix', required=True)
    parser.add_argument('--genomes', required=True)
    parser.add_argument('--fastani', default='fastANI')
    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.matrix, sep='\t', low_memory=False)
    print(f"Loaded {len(df)} pairs from {args.matrix}")

    tasks = [{
        'query': row['query'], 'reference': row['reference'],
        'genomes_dir': Path(args.genomes), 'fastani_path': args.fastani
    } for _, row in df.iterrows()]

    print(f"Running fastANI on {len(tasks)} pairs with {args.threads} threads...")
    with ProcessPoolExecutor(max_workers=args.threads) as pool:
        results = list(pool.map(run_fastani_pair, tasks))

    fastani_df = pd.DataFrame(results)
    merged = df.drop(columns=['fastani_ani'], errors='ignore').merge(
        fastani_df[['query', 'reference', 'fastani_ani']],
        on=['query', 'reference'], how='left'
    )

    # Reorder columns to put fastani_ani in original position if possible
    cols = list(df.columns)
    if 'fastani_ani' in cols:
        cols.remove('fastani_ani')
        idx = 11  # original position after r_species
        cols.insert(idx, 'fastani_ani')
        merged = merged[cols]

    merged.to_csv(args.output, sep='\t', index=False, float_format='%.6f')
    print(f"Saved {args.output}")
    print(f"fastANI missing: {merged['fastani_ani'].isna().sum()}/{len(merged)}")


if __name__ == '__main__':
    main()
