#!/usr/bin/env python3
"""Compare Syn2bANI performance across different enzymes and multi-enzyme mode.

Runs syn2bani dist --raw-features on a set of genome pairs using different
enzyme configurations, and records raw ANI, corrected ANI, shared tags,
alignment fractions, etc.
"""
import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd


SYN2BANI = Path('/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani')
GENOMES_DIR = Path('/lustre1/g/aos_shihuang/data/validation_mid_ani/genomes')


def run_syn2bani_pair(query, reference, enzyme=None, multi_enzyme=False):
    """Run syn2bani dist --raw-features on one pair."""
    q_path = GENOMES_DIR / f'{query}.fna'
    r_path = GENOMES_DIR / f'{reference}.fna'

    cmd = [str(SYN2BANI), 'dist', str(q_path), str(r_path),
           '--raw-features', '--min-af', '0.0']
    if multi_enzyme:
        cmd.append('--multi-enzyme')
    elif enzyme:
        cmd.extend(['-e', enzyme])

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            return {'error': result.stderr[:200]}

        # Parse raw-features TSV
        # query_file, ref_file, query_name, ref_name, raw_ani, af_q, af_r,
        # shared_tags, containment, div_proxy, ref_gc, corrected_ani
        for line in result.stdout.strip().split('\n'):
            if line.startswith('#') or not line.strip() or line.startswith('query_file'):
                continue
            parts = line.split('\t')
            if len(parts) >= 12:
                return {
                    'raw_ani': float(parts[4]),
                    'af_q': float(parts[5]),
                    'af_r': float(parts[6]),
                    'shared_tags': int(parts[7]),
                    'containment': float(parts[8]),
                    'ref_gc': float(parts[10]),
                    'corrected_ani': float(parts[11]),
                }
        return {'error': 'no_data_line'}
    except Exception as e:
        return {'error': str(e)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pairs', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    pairs = pd.read_csv(args.pairs, sep='\t')

    enzyme_modes = [
        ('BcgI', False),
        ('CjePI', False),
        ('CjeI', False),
        ('BsaXI', False),
        ('FalI', False),
        ('HaeIV', False),
        ('multi', True),
    ]

    records = []
    for _, row in pairs.iterrows():
        query = row['query']
        reference = row['reference']
        print(f'Processing {query} vs {reference}')

        for enzyme, multi in enzyme_modes:
            result = run_syn2bani_pair(query, reference, enzyme=enzyme if enzyme != 'multi' else None, multi_enzyme=multi)
            rec = {
                'query': query,
                'reference': reference,
                'q_species': row.get('q_species'),
                'r_species': row.get('r_species'),
                'q_genus': row.get('q_genus'),
                'r_genus': row.get('r_genus'),
                'enzyme_mode': enzyme,
            }
            if 'error' in result:
                rec['error'] = result['error']
            else:
                rec.update(result)
            records.append(rec)

    df = pd.DataFrame(records)
    df.to_csv(args.output, sep='\t', index=False, float_format='%.6f')
    print(f'Wrote {len(df)} records to {args.output}')


if __name__ == '__main__':
    main()
