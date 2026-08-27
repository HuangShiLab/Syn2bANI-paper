#!/usr/bin/env python3
"""Split a pair TSV into chunks for HPC job arrays."""
import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pairs', required=True, help='Input pair TSV')
    parser.add_argument('--output-dir', required=True, help='Directory for chunk files')
    parser.add_argument('--chunk-size', type=int, default=1000,
                        help='Pairs per chunk')
    parser.add_argument('--prefix', default='pairs_chunk',
                        help='Chunk file prefix')
    args = parser.parse_args()

    df = pd.read_csv(args.pairs, sep='\t')
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    n = len(df)
    n_chunks = (n + args.chunk_size - 1) // args.chunk_size
    for i in range(n_chunks):
        start = i * args.chunk_size
        end = min((i + 1) * args.chunk_size, n)
        chunk = df.iloc[start:end]
        chunk_path = out_dir / f'{args.prefix}_{i+1}.tsv'
        chunk.to_csv(chunk_path, sep='\t', index=False)

    print(f'Split {n} pairs into {n_chunks} chunks of <= {args.chunk_size} in {out_dir}')


if __name__ == '__main__':
    main()
