#!/usr/bin/env python3
"""
run_benchmark_matrix.py
Run skani, FastANI, and Syn2bANI (raw + GBRT) on a set of genome pairs.
Parallelized with multiprocessing for Mac Studio.

Usage:
  python3 run_benchmark_matrix.py \
    --pairs pairs_20k.tsv \
    --genomes ~/data/gtdb-r207/genomes/ \
    --syn2bani ~/Syn2bANI/target/release/syn2bani \
    --skani $(which skani) \
    --fastani $(which fastANI) \
    --output results/matrix.tsv \
    --threads 16
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile
from multiprocessing import Pool
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Runner helpers
# ---------------------------------------------------------------------------

def run_skani_pair(args: dict) -> dict:
    """Run skani dist on one pair. Returns {query, reference, ani, align_frac}."""
    q, r, skani_path = args['query'], args['reference'], args['skani_path']
    q_path = args['genomes_dir'] / f"{q}.fna"
    r_path = args['genomes_dir'] / f"{r}.fna"

    # If exact filename doesn't exist, try wildcard
    if not q_path.exists():
        matches = list(args['genomes_dir'].glob(f"*{q}*.fna"))
        if matches:
            q_path = matches[0]
    if not r_path.exists():
        matches = list(args['genomes_dir'].glob(f"*{r}*.fna"))
        if matches:
            r_path = matches[0]

    if not q_path.exists() or not r_path.exists():
        return {'query': q, 'reference': r, 'skani_ani': None,
                'skani_align_frac': None, 'error': 'missing_genome'}

    try:
        result = subprocess.run(
            [skani_path, 'dist', str(q_path), str(r_path)],
            capture_output=True, text=True, timeout=120
        )
        # Parse skani output: ANI\tAlign_fraction
        lines = result.stdout.strip().split('\n')
        for line in lines:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                return {
                    'query': q, 'reference': r,
                    'skani_ani': float(parts[0]),
                    'skani_align_frac': float(parts[1]) if len(parts) > 1 else None,
                    'error': None
                }
        return {'query': q, 'reference': r, 'skani_ani': None,
                'skani_align_frac': None, 'error': 'parse_failed'}
    except Exception as e:
        return {'query': q, 'reference': r, 'skani_ani': None,
                'skani_align_frac': None, 'error': str(e)}


def run_fastani_pair(args: dict) -> dict:
    """Run fastANI on one pair. Returns {query, reference, ani, frag_mapped, total_frag}."""
    q, r, fastani_path = args['query'], args['reference'], args['fastani_path']
    q_path = args['genomes_dir'] / f"{q}.fna"
    r_path = args['genomes_dir'] / f"{r}.fna"

    if not q_path.exists():
        matches = list(args['genomes_dir'].glob(f"*{q}*.fna"))
        if matches:
            q_path = matches[0]
    if not r_path.exists():
        matches = list(args['genomes_dir'].glob(f"*{r}*.fna"))
        if matches:
            r_path = matches[0]

    if not q_path.exists() or not r_path.exists():
        return {'query': q, 'reference': r, 'fastani_ani': None,
                'fastani_frag_mapped': None, 'fastani_total_frag': None,
                'error': 'missing_genome'}

    try:
        result = subprocess.run(
            [fastani_path, '-q', str(q_path), '-r', str(r_path),
             '-o', '/dev/stdout', '--noFrag'],
            capture_output=True, text=True, timeout=300
        )
        # fastANI output: q\tr\tANI\tfrag_mapped\ttotal_frag
        lines = result.stdout.strip().split('\n')
        for line in lines:
            if not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) >= 3:
                return {
                    'query': q, 'reference': r,
                    'fastani_ani': float(parts[2]),
                    'fastani_frag_mapped': int(parts[3]) if len(parts) > 3 else None,
                    'fastani_total_frag': int(parts[4]) if len(parts) > 4 else None,
                    'error': None
                }
        return {'query': q, 'reference': r, 'fastani_ani': None,
                'fastani_frag_mapped': None, 'fastani_total_frag': None,
                'error': 'parse_failed'}
    except Exception as e:
        return {'query': q, 'reference': r, 'fastani_ani': None,
                'fastani_frag_mapped': None, 'fastani_total_frag': None,
                'error': str(e)}


def run_syn2bani_pair(args: dict) -> dict:
    """Run Syn2bANI dist on one pair (raw + GBRT if available)."""
    q, r, syn2bani_path = args['query'], args['reference'], args['syn2bani_path']
    q_path = args['genomes_dir'] / f"{q}.fna"
    r_path = args['genomes_dir'] / f"{r}.fna"

    if not q_path.exists():
        matches = list(args['genomes_dir'].glob(f"*{q}*.fna"))
        if matches:
            q_path = matches[0]
    if not r_path.exists():
        matches = list(args['genomes_dir'].glob(f"*{r}*.fna"))
        if matches:
            r_path = matches[0]

    if not q_path.exists() or not r_path.exists():
        return {'query': q, 'reference': r, 's2b_raw_ani': None,
                's2b_gbrt_ani': None, 's2b_shared_tags': None,
                'error': 'missing_genome'}

    results = {'query': q, 'reference': r}

    # --- Raw Syn2bANI ---
    try:
        result = subprocess.run(
            [syn2bani_path, 'dist', str(q_path), str(r_path),
             '-e', 'BcgI', '--no-debias'],
            capture_output=True, text=True, timeout=120
        )
        # Parse TSV output
        for line in result.stdout.strip().split('\n'):
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) >= 6:
                results['s2b_raw_ani'] = float(parts[5])
                results['s2b_shared_tags'] = int(parts[6]) if len(parts) > 6 else None
                break
    except Exception as e:
        results['s2b_raw_ani'] = None
        results['s2b_shared_tags'] = None
        results['error_raw'] = str(e)

    # --- GBRT debiased Syn2bANI ---
    try:
        result = subprocess.run(
            [syn2bani_path, 'dist', str(q_path), str(r_path),
             '-e', 'BcgI'],
            capture_output=True, text=True, timeout=120
        )
        for line in result.stdout.strip().split('\n'):
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) >= 6:
                results['s2b_gbrt_ani'] = float(parts[5])
                break
    except Exception as e:
        results['s2b_gbrt_ani'] = None
        results['error_gbrt'] = str(e)

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pairs', required=True)
    parser.add_argument('--genomes', required=True)
    parser.add_argument('--syn2bani', default='syn2bani')
    parser.add_argument('--skani', default='skani')
    parser.add_argument('--fastani', default='fastANI')
    parser.add_argument('--output', required=True)
    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--tools', default='all',
                        help='Comma-separated: skani,fastani,syn2bani,all')
    args = parser.parse_args()

    genomes_dir = Path(args.genomes)
    if not genomes_dir.exists():
        print(f"ERROR: Genomes directory not found: {genomes_dir}")
        sys.exit(1)

    pairs_df = pd.read_csv(args.pairs, sep='\t')
    print(f"Loaded {len(pairs_df)} pairs from {args.pairs}")

    tools = [t.strip() for t in args.tools.split(',')]
    run_all = 'all' in tools

    # --- Run skani ---
    if run_all or 'skani' in tools:
        print(f"\n[1/3] Running skani on {len(pairs_df)} pairs...")
        skani_args = [{
            'query': row['query'], 'reference': row['reference'],
            'genomes_dir': genomes_dir, 'skani_path': args.skani
        } for _, row in pairs_df.iterrows()]

        with Pool(args.threads) as pool:
            skani_results = pool.map(run_skani_pair, skani_args)

        skani_df = pd.DataFrame(skani_results)
        pairs_df = pairs_df.merge(
            skani_df[['query', 'reference', 'skani_ani', 'skani_align_frac']],
            on=['query', 'reference'], how='left'
        )
        print(f"  skani completed. Missing: {pairs_df['skani_ani'].isna().sum()}")

    # --- Run FastANI ---
    if run_all or 'fastani' in tools:
        print(f"\n[2/3] Running FastANI on {len(pairs_df)} pairs...")
        fastani_args = [{
            'query': row['query'], 'reference': row['reference'],
            'genomes_dir': genomes_dir, 'fastani_path': args.fastani
        } for _, row in pairs_df.iterrows()]

        with Pool(args.threads) as pool:
            fastani_results = pool.map(run_fastani_pair, fastani_args)

        fastani_df = pd.DataFrame(fastani_results)
        pairs_df = pairs_df.merge(
            fastani_df[['query', 'reference', 'fastani_ani',
                        'fastani_frag_mapped', 'fastani_total_frag']],
            on=['query', 'reference'], how='left'
        )
        print(f"  FastANI completed. Missing: {pairs_df['fastani_ani'].isna().sum()}")

    # --- Run Syn2bANI ---
    if run_all or 'syn2bani' in tools:
        print(f"\n[3/3] Running Syn2bANI on {len(pairs_df)} pairs...")
        s2b_args = [{
            'query': row['query'], 'reference': row['reference'],
            'genomes_dir': genomes_dir, 'syn2bani_path': args.syn2bani
        } for _, row in pairs_df.iterrows()]

        with Pool(args.threads) as pool:
            s2b_results = pool.map(run_syn2bani_pair, s2b_args)

        s2b_df = pd.DataFrame(s2b_results)
        pairs_df = pairs_df.merge(
            s2b_df[['query', 'reference', 's2b_raw_ani',
                    's2b_gbrt_ani', 's2b_shared_tags']],
            on=['query', 'reference'], how='left'
        )
        print(f"  Syn2bANI completed. Missing raw: {pairs_df['s2b_raw_ani'].isna().sum()}")

    # --- Save ---
    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    pairs_df.to_csv(args.output, sep='\t', index=False, float_format='%.4f')
    print(f"\n{'='*60}")
    print(f"Results saved: {args.output}")
    print(f"Columns: {list(pairs_df.columns)}")
    print(f"{'='*60}")
    print("\nNext steps:")
    print(f"  python3 scripts/train_gbrt_v3.py --matrix {args.output}")
    print(f"  python3 scripts/evaluate_by_taxonomy.py --matrix {args.output}")


if __name__ == '__main__':
    main()
