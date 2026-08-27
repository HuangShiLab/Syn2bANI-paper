#!/usr/bin/env python3
"""
Time each tool on the 1k GTDB-R207 pair sample.
Uses the same per-pair runners as run_benchmark_matrix_v2.py.
"""
import time
import sys
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

# Re-use runners from the benchmark matrix script
sys.path.insert(0, str(Path(__file__).parent))
from run_benchmark_matrix_v2 import (
    resolve_genome_path,
    run_fastani_pair_pyfastani,
    run_skani_pair,
    run_syn2bani_pair,
)

PAIRS = Path('results/pairs_gtdb_r207_1k.tsv')
GENOMES = Path.home() / 'data' / 'gtdb-r207' / 'genomes_all'
SYN2BANI = Path.home() / 'Downloads' / 'Syn2bANI' / 'target' / 'release' / 'syn2bani'
SKANI = 'skani'
THREADS = 16


def time_phase(name, runner, args, n_pairs):
    print(f"\nTiming {name} on {n_pairs} pairs...")
    t0 = time.perf_counter()
    results = [runner(a) for a in args]
    elapsed = time.perf_counter() - t0
    print(f"  {name}: {elapsed:.2f}s ({n_pairs/elapsed:.1f} pairs/s)")
    missing = sum(1 for r in results if r.get('fastani_ani') is None
                  and r.get('skani_ani') is None
                  and r.get('s2b_raw_ani') is None)
    print(f"  Missing results: {missing}")
    return elapsed


def time_phase_parallel(name, runner, args, n_pairs, threads):
    print(f"\nTiming {name} on {n_pairs} pairs ({threads} processes)...")
    t0 = time.perf_counter()
    with Pool(threads) as pool:
        results = pool.map(runner, args)
    elapsed = time.perf_counter() - t0
    print(f"  {name}: {elapsed:.2f}s ({n_pairs/elapsed:.1f} pairs/s)")
    missing = sum(1 for r in results if r.get('skani_ani') is None
                  and r.get('s2b_raw_ani') is None)
    print(f"  Missing results: {missing}")
    return elapsed


def main():
    pairs_df = pd.read_csv(PAIRS, sep='\t')
    n = len(pairs_df)
    print(f"Pairs: {n}")

    fastani_args = [{
        'query': row['query'], 'reference': row['reference'],
        'genomes_dir': GENOMES, 'fastani_path': 'fastANI'
    } for _, row in pairs_df.iterrows()]

    skani_args = [{
        'query': row['query'], 'reference': row['reference'],
        'genomes_dir': GENOMES, 'skani_path': SKANI
    } for _, row in pairs_df.iterrows()]

    s2b_args = [{
        'query': row['query'], 'reference': row['reference'],
        'genomes_dir': GENOMES, 'syn2bani_path': str(SYN2BANI)
    } for _, row in pairs_df.iterrows()]

    # FastANI is sequential in the matrix runner
    t_fastani = time_phase('FastANI (pyfastani, sequential)',
                           run_fastani_pair_pyfastani, fastani_args, n)

    # skani parallel
    t_skani = time_phase_parallel('skani', run_skani_pair, skani_args, n, THREADS)

    # Syn2bANI parallel
    t_s2b = time_phase_parallel('Syn2bANI optimized', run_syn2bani_pair, s2b_args, n, THREADS)

    print("\n" + "=" * 60)
    print("Summary (1,000 pairs, wall time)")
    print(f"  FastANI (pyfastani, sequential): {t_fastani:.2f}s")
    print(f"  skani ({THREADS} proc):            {t_skani:.2f}s")
    print(f"  Syn2bANI optimized ({THREADS} proc): {t_s2b:.2f}s")
    print("=" * 60)


if __name__ == '__main__':
    main()
