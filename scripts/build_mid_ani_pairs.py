#!/usr/bin/env python3
"""Build within-genus, different-species pairs for mid-ANI validation.

Only pairs from the same genus but different species are retained, because
this is the ANI range (approximately 85-95%) where Syn2bANI GBRT v4 needs
independent validation against FastANI.
"""
import argparse
import itertools

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description='Build within-genus, different-species validation pairs.'
    )
    parser.add_argument('--manifest', required=True,
                        help='TSV manifest with accession, species, genus, file columns')
    parser.add_argument('--output', required=True,
                        help='Output TSV of pairs')
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep='\t')
    required = {'accession', 'species', 'genus'}
    missing = required - set(manifest.columns)
    if missing:
        raise ValueError(f'Manifest missing columns: {missing}')

    records = manifest.to_dict('records')

    pairs = []
    for a, b in itertools.combinations(records, 2):
        if a['genus'] != b['genus']:
            continue
        if a['species'] == b['species']:
            continue
        pairs.append({
            'query': a['accession'],
            'reference': b['accession'],
            'label': 'mid_high',
            'q_species': a['species'],
            'r_species': b['species'],
            'q_genus': a['genus'],
            'r_genus': b['genus'],
            'q_category': 'mid_ani_validation',
            'r_category': 'mid_ani_validation',
        })

    df = pd.DataFrame(pairs)
    df.to_csv(args.output, sep='\t', index=False)
    print(f'Built {len(df)} within-genus, different-species pairs: {args.output}')
    if len(df) > 0:
        print(df['q_genus'].value_counts())


if __name__ == '__main__':
    main()
