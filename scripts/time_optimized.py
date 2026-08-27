#!/usr/bin/env python3
"""Time the optimized Syn2bANI binary on GTDB-R207 pairs."""
import time
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
import sys
sys.path.insert(0, str(Path(__file__).parent))
from run_benchmark_matrix_v2 import resolve_genome_path, run_syn2bani_pair


def time_s2b(pairs_tsv, binary, threads=16):
    pairs_df = pd.read_csv(pairs_tsv, sep='\t')
    n = len(pairs_df)
    print(f'Timing {binary} on {n} pairs...')
    genomes = Path.home() / 'data' / 'gtdb-r207' / 'genomes_all'
    args = [{
        'query': row['query'], 'reference': row['reference'],
        'genomes_dir': genomes, 'syn2bani_path': binary
    } for _, row in pairs_df.iterrows()]
    t0 = time.perf_counter()
    with Pool(threads) as pool:
        pool.map(run_syn2bani_pair, args)
    elapsed = time.perf_counter() - t0
    print(f'  {binary} {n} pairs: {elapsed:.2f}s ({n/elapsed:.1f} pairs/s)')
    return elapsed


if __name__ == '__main__':
    opt = Path.home() / 'Downloads' / 'Syn2bANI' / 'target' / 'release' / 'syn2bani'
    time_s2b('results/pairs_gtdb_r207_100.tsv', str(opt), threads=16)
    time_s2b('results/pairs_gtdb_r207_1k.tsv', str(opt), threads=16)
