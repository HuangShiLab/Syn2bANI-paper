#!/usr/bin/env python3
"""Time skani and FastANI on the 100-pair GTDB-R207 sample."""
import time
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
import sys
sys.path.insert(0, str(Path(__file__).parent))
from run_benchmark_matrix_v2 import run_skani_pair, run_fastani_pair_pyfastani


def time_skani(pairs_tsv, threads=16):
    pairs_df = pd.read_csv(pairs_tsv, sep='\t')
    n = len(pairs_df)
    genomes = Path.home() / 'data' / 'gtdb-r207' / 'genomes_all'
    args = [{'query': row['query'], 'reference': row['reference'],
             'genomes_dir': genomes, 'skani_path': 'skani'} for _, row in pairs_df.iterrows()]
    t0 = time.perf_counter()
    with Pool(threads) as pool:
        pool.map(run_skani_pair, args)
    elapsed = time.perf_counter() - t0
    print(f'skani {n} pairs: {elapsed:.2f}s ({n/elapsed:.1f} pairs/s)')


def time_fastani(pairs_tsv):
    pairs_df = pd.read_csv(pairs_tsv, sep='\t')
    n = len(pairs_df)
    genomes = Path.home() / 'data' / 'gtdb-r207' / 'genomes_all'
    args = [{'query': row['query'], 'reference': row['reference'],
             'genomes_dir': genomes, 'fastani_path': 'fastANI'} for _, row in pairs_df.iterrows()]
    t0 = time.perf_counter()
    for a in args:
        run_fastani_pair_pyfastani(a)
    elapsed = time.perf_counter() - t0
    print(f'FastANI (pyfastani, sequential) {n} pairs: {elapsed:.2f}s ({n/elapsed:.1f} pairs/s)')


if __name__ == '__main__':
    time_skani('results/pairs_gtdb_r207_100.tsv', threads=16)
    time_fastani('results/pairs_gtdb_r207_100.tsv')
