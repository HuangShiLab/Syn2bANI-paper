#!/usr/bin/env python3
"""Prepare validation matrix with BcgI, CjePI, or combined features."""
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bcgi-matrix', required=True)
    parser.add_argument('--cjepi-matrix', required=True)
    parser.add_argument('--reference-matrix', required=True,
                        help='Matrix with fastani_ani reference')
    parser.add_argument('--mode', required=True, choices=['bcgi', 'cjepi', 'combined'])
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    bcgi = pd.read_csv(args.bcgi_matrix, sep='\t', low_memory=False)
    cjepi = pd.read_csv(args.cjepi_matrix, sep='\t', low_memory=False)
    ref = pd.read_csv(args.reference_matrix, sep='\t', low_memory=False)[['query', 'reference', 'fastani_ani', 'skani_ani']]

    # BcgI matrix: rename relevant columns to bcgi_ prefix
    bcgi_cols = {
        's2b_raw_ani': 'bcgi_raw_ani',
        's2b_mash_ani': 'bcgi_mash_ani',
        's2b_af_q': 'bcgi_af_q',
        's2b_af_r': 'bcgi_af_r',
        's2b_shared_tags': 'bcgi_shared_tags',
        's2b_gbrt_ani': 'bcgi_gbrt_ani',
    }
    bcgi = bcgi.rename(columns=bcgi_cols)

    # CjePI matrix already has cjepi_ prefix
    if 'error' in cjepi.columns:
        cjepi = cjepi[cjepi['error'].isna()]

    # Merge
    merged = bcgi.merge(cjepi, on=['query', 'reference'], how='inner' if args.mode == 'combined' else 'left')
    merged = merged.merge(ref, on=['query', 'reference'], how='left')

    keep = ['query', 'reference', 'fastani_ani', 'skani_ani', 'label']
    if 'q_species' in merged.columns:
        keep.extend(['q_species', 'r_species', 'q_genus', 'r_genus'])

    if args.mode in ('bcgi', 'combined'):
        keep.extend(['bcgi_raw_ani', 'bcgi_mash_ani', 'bcgi_af_q', 'bcgi_af_r', 'bcgi_shared_tags'])
    if args.mode in ('cjepi', 'combined'):
        keep.extend(['cjepi_raw_ani', 'cjepi_mash_ani', 'cjepi_af_q', 'cjepi_af_r', 'cjepi_shared_tags'])

    merged = merged[[c for c in keep if c in merged.columns]]

    merged.to_csv(args.output, sep='\t', index=False, float_format='%.6f')
    print(f'Wrote {len(merged)} rows to {args.output} (mode={args.mode})')


if __name__ == '__main__':
    main()
