#!/usr/bin/env python3
"""Prepare training matrix with BcgI, CjePI, or combined features."""
import argparse
import numpy as np
import pandas as pd


def compute_mash_ani(af_q, af_r, tag_len=32.0):
    c1 = np.maximum(af_q, 1e-10)
    c2 = np.maximum(af_r, 1e-10)
    containment_geo = np.sqrt(c1 * c2)
    return np.clip(1.0 + np.log(containment_geo) / tag_len, 0.0, 1.0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bcgi-matrix', required=True)
    parser.add_argument('--cjepi-matrix', required=True)
    parser.add_argument('--mode', required=True, choices=['bcgi', 'cjepi', 'combined'])
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    bcgi = pd.read_csv(args.bcgi_matrix, sep='\t', low_memory=False)
    cjepi = pd.read_csv(args.cjepi_matrix, sep='\t', low_memory=False)

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
    if 'bcgi_mash_ani' not in bcgi.columns:
        bcgi['bcgi_mash_ani'] = compute_mash_ani(bcgi['bcgi_af_q'].values, bcgi['bcgi_af_r'].values)

    # CjePI matrix already has cjepi_ prefix from run_benchmark_enzyme.py
    # But it has error column possibly; keep only successful rows
    if 'error' in cjepi.columns:
        cjepi = cjepi[cjepi['error'].isna()]
    if 'cjepi_mash_ani' not in cjepi.columns:
        cjepi['cjepi_mash_ani'] = compute_mash_ani(cjepi['cjepi_af_q'].values, cjepi['cjepi_af_r'].values)

    # Merge on query/reference
    merged = bcgi.merge(cjepi, on=['query', 'reference'], how='inner' if args.mode == 'combined' else 'left')

    # Keep all relevant columns
    keep = ['query', 'reference', 'fastani_ani', 'label']
    if 'q_species' in merged.columns:
        keep.extend(['q_species', 'r_species', 'q_genus', 'r_genus'])

    if args.mode in ('bcgi', 'combined'):
        keep.extend(['bcgi_raw_ani', 'bcgi_mash_ani', 'bcgi_af_q', 'bcgi_af_r', 'bcgi_shared_tags'])
    if args.mode in ('cjepi', 'combined'):
        keep.extend(['cjepi_raw_ani', 'cjepi_mash_ani', 'cjepi_af_q', 'cjepi_af_r', 'cjepi_shared_tags'])

    merged = merged[[c for c in keep if c in merged.columns]]
    merged = merged.dropna(subset=['fastani_ani'])

    merged.to_csv(args.output, sep='\t', index=False, float_format='%.6f')
    print(f'Wrote {len(merged)} rows to {args.output} (mode={args.mode})')


if __name__ == '__main__':
    main()
