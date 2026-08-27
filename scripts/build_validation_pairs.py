#!/usr/bin/env python3
"""Build all-vs-all validation pairs from oral/gut genomes."""
import argparse
import itertools
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep='\t')
    records = manifest.to_dict('records')

    # Extract genus from species name (first word)
    for r in records:
        r['genus'] = r['species'].split()[0]
        r['phylum'] = 'unknown'  # will be filled later if needed

    pairs = []
    for a, b in itertools.combinations(records, 2):
        # Use accession as identifier
        q = a['accession']
        r = b['accession']

        # Determine label based on taxonomy
        if a['species'] == b['species']:
            label = 'high'
        elif a['genus'] == b['genus']:
            label = 'mid_high'
        else:
            label = 'low'

        pairs.append({
            'query': q,
            'reference': r,
            'label': label,
            'q_species': a['species'],
            'r_species': b['species'],
            'q_genus': a['genus'],
            'r_genus': b['genus'],
            'q_category': a['category'],
            'r_category': b['category'],
        })

    df = pd.DataFrame(pairs)
    df.to_csv(args.output, sep='\t', index=False)
    print(f'Built {len(df)} pairs: {args.output}')
    print(df['label'].value_counts())


if __name__ == '__main__':
    main()
