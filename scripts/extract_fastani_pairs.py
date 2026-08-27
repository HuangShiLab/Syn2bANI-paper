#!/usr/bin/env python3
"""Extract pairs with FastANI reference from a matrix."""
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.matrix, sep='\t', low_memory=False)
    valid = df[df['fastani_ani'].notna()].copy()

    # Keep minimal columns needed for pair running
    keep = ['query', 'reference', 'label']
    if 'q_species' in valid.columns:
        keep.extend(['q_species', 'r_species', 'q_genus', 'r_genus'])
    valid = valid[[c for c in keep if c in valid.columns]]

    valid.to_csv(args.output, sep='\t', index=False)
    print(f'Extracted {len(valid)} pairs with FastANI reference to {args.output}')


if __name__ == '__main__':
    main()
