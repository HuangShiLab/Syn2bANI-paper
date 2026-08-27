#!/usr/bin/env python3
"""Run Syn2bANI with all 16 single enzymes and multi-enzyme mode on validation pairs.

Parses the raw-features TSV correctly, including the mash_ani column.
"""
import argparse
import subprocess
from pathlib import Path

import pandas as pd


SYN2BANI = Path('/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani')
GENOMES_DIR = Path('/lustre1/g/aos_shihuang/data/validation_mid_ani/genomes')
TIMEOUT = 600  # seconds

ENZYMES = [
    'BcgI', 'AlfI', 'AloI', 'BaeI', 'BplI', 'BsaXI', 'BslFI', 'Bsp24I',
    'CjeI', 'CjePI', 'CspCI', 'FalI', 'HaeIV', 'Hin4I', 'PpiI', 'PsrI',
]


def run_pair(query, reference, enzyme=None, multi_enzyme=False):
    q_path = GENOMES_DIR / f'{query}.fna'
    r_path = GENOMES_DIR / f'{reference}.fna'
    cmd = [str(SYN2BANI), 'dist', str(q_path), str(r_path),
           '--raw-features', '--min-af', '0.0']
    if multi_enzyme:
        cmd.append('--multi-enzyme')
    elif enzyme:
        cmd.extend(['-e', enzyme])

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT)
        if result.returncode != 0:
            return {'error': result.stderr[:300]}
        for line in result.stdout.strip().split('\n'):
            if line.startswith('#') or not line.strip() or line.startswith('query_file'):
                continue
            parts = line.split('\t')
            if len(parts) >= 13:
                return {
                    'raw_ani': float(parts[4]),
                    'mash_ani': float(parts[5]),
                    'af_q': float(parts[6]),
                    'af_r': float(parts[7]),
                    'shared_tags': int(parts[8]),
                    'containment': float(parts[9]),
                    'ref_gc': float(parts[11]),
                    'corrected_ani': float(parts[12]),
                }
        return {'error': 'no_data_line'}
    except subprocess.TimeoutExpired:
        return {'error': f'timeout_after_{TIMEOUT}s'}
    except Exception as e:
        return {'error': str(e)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pairs', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    pairs = pd.read_csv(args.pairs, sep='\t')
    records = []

    modes = [('multi', None, True)]
    for enz in ENZYMES:
        modes.append((enz, enz, False))

    total = len(pairs) * len(modes)
    done = 0
    for _, row in pairs.iterrows():
        query = row['query']
        reference = row['reference']
        for mode_name, enzyme, multi in modes:
            done += 1
            print(f'[{done}/{total}] {query} vs {reference} :: {mode_name}', flush=True)
            result = run_pair(query, reference, enzyme=enzyme, multi_enzyme=multi)
            rec = {
                'query': query,
                'reference': reference,
                'q_species': row.get('q_species'),
                'r_species': row.get('r_species'),
                'q_genus': row.get('q_genus'),
                'r_genus': row.get('r_genus'),
                'enzyme_mode': mode_name,
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
