#!/usr/bin/env python3
"""Generate multi-enzyme features for training pairs."""
import argparse
import subprocess
from multiprocessing import Pool
from pathlib import Path

import pandas as pd


SYN2BANI = Path('/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani')
GENOMES = Path('/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all')
TIMEOUT = 120


def run_pair(args):
    idx, query, reference = args
    q_path = GENOMES / f'{query}.fna'
    r_path = GENOMES / f'{reference}.fna'
    cmd = [str(SYN2BANI), 'dist', str(q_path), str(r_path),
           '--raw-features', '--min-af', '0.0', '--multi-enzyme']
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT)
        if res.returncode != 0:
            return {'query': query, 'reference': reference, 'error': res.stderr[:200]}
        for line in res.stdout.strip().split('\n'):
            if line.startswith('#') or not line.strip() or line.startswith('query_file'):
                continue
            parts = line.split('\t')
            if len(parts) >= 13:
                return {
                    'query': query,
                    'reference': reference,
                    'multi_raw_ani': float(parts[4]),
                    'multi_mash_ani': float(parts[5]),
                    'multi_af_q': float(parts[6]),
                    'multi_af_r': float(parts[7]),
                    'multi_shared_tags': int(parts[8]),
                }
        return {'query': query, 'reference': reference, 'error': 'no_data'}
    except subprocess.TimeoutExpired:
        return {'query': query, 'reference': reference, 'error': 'timeout'}
    except Exception as e:
        return {'query': query, 'reference': reference, 'error': str(e)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pairs', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--threads', type=int, default=16)
    args = parser.parse_args()

    pairs = pd.read_csv(args.pairs, sep='\t')
    tasks = [(i, row['query'], row['reference']) for i, row in pairs.iterrows()]

    records = []
    with Pool(processes=args.threads) as pool:
        for rec in pool.imap_unordered(run_pair, tasks):
            records.append(rec)

    df = pd.DataFrame(records)
    df = pairs.merge(df, on=['query', 'reference'], how='left')
    df.to_csv(args.output, sep='\t', index=False, float_format='%.6f')
    print(f'Wrote {len(df)} records to {args.output}')
    if 'error' in df.columns:
        n_err = df['error'].notna().sum()
        print(f'Errors: {n_err}')


if __name__ == '__main__':
    main()
