#!/usr/bin/env python3
"""
run_benchmark_matrix_v8.py
Run skani and Syn2bANI v8 (ani path) on a set of genome pairs.

Supports checkpoint/resume: if --output exists, already-completed pairs are skipped.

Usage:
  python3 scripts/run_benchmark_matrix_v8.py \
    --pairs results/pairs_gtdb_r207.tsv \
    --genomes ~/data/gtdb-r207/genomes_all/ \
    --syn2bani ~/Downloads/Syn2bANI/target/release/syn2bani \
    --skani $(which skani) \
    --output results/matrix.tsv \
    --threads 16
"""

import argparse
import os
import subprocess
import sys
import time
from multiprocessing import Pool
from pathlib import Path
from typing import Optional

import pandas as pd
import numpy as np


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def resolve_genome_path(accession: str, genomes_dir) -> Optional[Path]:
    """Resolve an accession to a genome file in genomes_dir."""
    genomes_dir = Path(genomes_dir)
    direct = genomes_dir / f"{accession}.fna"
    if direct.exists():
        return direct
    matches = list(genomes_dir.glob(f"*{accession}*.fna"))
    if matches:
        return matches[0]
    return None


# ---------------------------------------------------------------------------
# Syn2bANI v8 ani
# ---------------------------------------------------------------------------

S2B_V8_COLUMNS = [
    'query', 'reference', 'ani', 'ani_uniform', 'af_query', 'af_reference',
    'std_err', 'het_shape', 'retention', 'ani_from_loss', 'ani_from_hist',
    'n_anchors', 'n_chains', 'n_tags', 'flag'
]


def _parse_syn2bani_v8_line(stdout: str) -> Optional[dict]:
    """Parse the single data line from `syn2bani ani --verbose`."""
    for line in stdout.strip().split('\n'):
        if line.startswith('query\treference'):
            continue
        if not line.strip():
            continue
        parts = line.split('\t')
        if len(parts) < 15:
            continue
        try:
            return {
                's2b_ani': float(parts[2]),
                's2b_ani_uniform': float(parts[3]),
                's2b_af_q': float(parts[4]),
                's2b_af_r': float(parts[5]),
                's2b_std_err': float(parts[6]),
                's2b_het_shape': parts[7],
                's2b_retention': float(parts[8]),
                's2b_ani_from_loss': float(parts[9]),
                's2b_ani_from_hist': float(parts[10]),
                's2b_n_anchors': int(parts[11]),
                's2b_n_chains': int(parts[12]),
                's2b_n_tags': int(parts[13]),
                's2b_flag': parts[14],
            }
        except (ValueError, IndexError):
            return None
    return None


def run_syn2bani_v8_pair(args: dict) -> dict:
    """Run Syn2bANI v8 ani on one pair."""
    q, r, syn2bani_path = args['query'], args['reference'], args['syn2bani_path']
    q_path = resolve_genome_path(q, args['genomes_dir'])
    r_path = resolve_genome_path(r, args['genomes_dir'])

    results = {'query': q, 'reference': r}
    if q_path is None or r_path is None:
        results['error'] = 'missing_genome'
        for col in ['s2b_ani', 's2b_ani_uniform', 's2b_af_q', 's2b_af_r',
                    's2b_std_err', 's2b_het_shape', 's2b_retention',
                    's2b_ani_from_loss', 's2b_ani_from_hist',
                    's2b_n_anchors', 's2b_n_chains', 's2b_n_tags', 's2b_flag']:
            results[col] = None
        return results

    try:
        cmd = [syn2bani_path, 'ani', str(q_path), str(r_path), '--verbose', '-t', '1']
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=120
        )
        parsed = _parse_syn2bani_v8_line(result.stdout)
        if parsed:
            results.update(parsed)
            results['error'] = None
        else:
            for col in ['s2b_ani', 's2b_ani_uniform', 's2b_af_q', 's2b_af_r',
                        's2b_std_err', 's2b_het_shape', 's2b_retention',
                        's2b_ani_from_loss', 's2b_ani_from_hist',
                        's2b_n_anchors', 's2b_n_chains', 's2b_n_tags', 's2b_flag']:
                results[col] = None
            results['error'] = f'parse_failed rc={result.returncode} stdout={result.stdout[:200]} stderr={result.stderr[:200]}'
    except Exception as e:
        for col in ['s2b_ani', 's2b_ani_uniform', 's2b_af_q', 's2b_af_r',
                    's2b_std_err', 's2b_het_shape', 's2b_retention',
                    's2b_ani_from_loss', 's2b_ani_from_hist',
                    's2b_n_anchors', 's2b_n_chains', 's2b_n_tags', 's2b_flag']:
            results[col] = None
        results['error'] = str(e)

    return results


# ---------------------------------------------------------------------------
# skani
# ---------------------------------------------------------------------------

def run_skani_pair(args: dict) -> dict:
    """Run skani dist on one pair."""
    q, r, skani_path = args['query'], args['reference'], args['skani_path']
    q_path = resolve_genome_path(q, args['genomes_dir'])
    r_path = resolve_genome_path(r, args['genomes_dir'])

    results = {'query': q, 'reference': r}
    if q_path is None or r_path is None:
        results['skani_ani'] = None
        results['skani_align_frac_ref'] = None
        results['skani_align_frac_query'] = None
        results['error'] = 'missing_genome'
        return results

    try:
        cmd = [skani_path, 'dist', '-q', str(q_path), '-r', str(r_path), '-t', '1']
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=120
        )
        for line in result.stdout.strip().split('\n'):
            if line.startswith('Ref_file') or not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) >= 5:
                try:
                    results['skani_ani'] = float(parts[2]) / 100.0
                    results['skani_align_frac_ref'] = float(parts[3]) / 100.0
                    results['skani_align_frac_query'] = float(parts[4]) / 100.0
                    results['error'] = None
                    return results
                except ValueError:
                    pass
        results['skani_ani'] = None
        results['skani_align_frac_ref'] = None
        results['skani_align_frac_query'] = None
        results['error'] = f'parse_failed rc={result.returncode} stdout={result.stdout[:200]}'
    except Exception as e:
        results['skani_ani'] = None
        results['skani_align_frac_ref'] = None
        results['skani_align_frac_query'] = None
        results['error'] = str(e)

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
    parser.add_argument('--output', required=True)
    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--tools', default='all',
                        help='Comma-separated: skani,syn2bani,all')
    parser.add_argument('--chunk-size', type=int, default=1000,
                        help='Save checkpoint every N pairs')
    args = parser.parse_args()

    genomes_dir = Path(args.genomes)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Load pairs
    pairs_df = pd.read_csv(args.pairs, sep='\t', low_memory=False)
    print(f"Loaded {len(pairs_df)} pairs from {args.pairs}")

    # Resume if output exists
    if output_path.exists():
        existing = pd.read_csv(output_path, sep='\t', low_memory=False)
        print(f"Resuming from {output_path}: {len(existing)} rows")
        # Determine which pairs are already done based on tools present
        done = existing.copy()
        if 'syn2bani' in args.tools or args.tools == 'all':
            done = done[done['s2b_ani'].notna()]
        if 'skani' in args.tools or args.tools == 'all':
            done = done[done['skani_ani'].notna()]
        done_keys = set(zip(done['query'], done['reference']))
        pairs_df = pairs_df[~pairs_df.apply(lambda row: (row['query'], row['reference']) in done_keys, axis=1)]
        print(f"Remaining pairs: {len(pairs_df)}")
    else:
        existing = pd.DataFrame()

    run_all = args.tools == 'all'
    tools = set(args.tools.split(',')) if args.tools != 'all' else {'skani', 'syn2bani'}

    remaining = pairs_df.copy()

    # --- Run skani ---
    if run_all or 'skani' in tools:
        print(f"\n[1/2] Running skani on {len(remaining)} pairs...")
        skani_args = [{
            'query': row['query'], 'reference': row['reference'],
            'genomes_dir': genomes_dir, 'skani_path': args.skani
        } for _, row in remaining.iterrows()]

        with Pool(args.threads) as pool:
            skani_results = pool.map(run_skani_pair, skani_args)

        skani_df = pd.DataFrame(skani_results)
        remaining = remaining.merge(
            skani_df[['query', 'reference', 'skani_ani', 'skani_align_frac_ref', 'skani_align_frac_query']],
            on=['query', 'reference'], how='left'
        )
        print(f"  skani completed. Missing: {remaining['skani_ani'].isna().sum()}")

    # --- Run Syn2bANI v8 ---
    if run_all or 'syn2bani' in tools:
        print(f"\n[2/2] Running Syn2bANI v8 on {len(remaining)} pairs...")
        s2b_args = [{
            'query': row['query'], 'reference': row['reference'],
            'genomes_dir': genomes_dir, 'syn2bani_path': args.syn2bani
        } for _, row in remaining.iterrows()]

        with Pool(args.threads) as pool:
            s2b_results = pool.map(run_syn2bani_v8_pair, s2b_args)

        s2b_df = pd.DataFrame(s2b_results)
        remaining = remaining.merge(
            s2b_df[['query', 'reference', 's2b_ani', 's2b_ani_uniform', 's2b_af_q', 's2b_af_r',
                    's2b_std_err', 's2b_het_shape', 's2b_retention',
                    's2b_ani_from_loss', 's2b_ani_from_hist',
                    's2b_n_anchors', 's2b_n_chains', 's2b_n_tags', 's2b_flag']],
            on=['query', 'reference'], how='left'
        )
        print(f"  Syn2bANI v8 completed. Missing: {remaining['s2b_ani'].isna().sum()}")

    # --- Merge with existing and save ---
    os.makedirs(output_path.parent or '.', exist_ok=True)
    if output_path.exists() and not existing.empty:
        combined = pd.concat([existing, remaining], ignore_index=True)
        combined = combined.drop_duplicates(subset=['query', 'reference'], keep='last')
    else:
        combined = remaining

    combined.to_csv(output_path, sep='\t', index=False, float_format='%.6f')
    print(f"\nSaved {len(combined)} rows to {output_path}")


if __name__ == '__main__':
    main()
