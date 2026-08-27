#!/usr/bin/env python3
"""Run Syn2bANI with all 16 single enzymes on validation pairs (parallel)."""
import argparse
import subprocess
from multiprocessing import Pool, cpu_count
from pathlib import Path

import pandas as pd


SYN2BANI = Path('/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani')
GENOMES_DIR = Path('/lustre1/g/aos_shihuang/data/validation_mid_ani/genomes')
TIMEOUT = 120

ENZYMES = [
    'BcgI', 'AlfI', 'AloI', 'BaeI', 'BplI', 'BsaXI', 'BslFI', 'Bsp24I',
    'CjeI', 'CjePI', 'CspCI', 'FalI', 'HaeIV', 'Hin4I', 'PpiI', 'PsrI',
]


def run_pair(query, reference, enzyme):
    q_path = GENOMES_DIR / f'{query}.fna'
    r_path = GENOMES_DIR / f'{reference}.fna'
    cmd = [str(SYN2BANI), 'dist', str(q_path), str(r_path),
           '--raw-features', '--min-af', '0.0', '-e', enzyme]
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


def process_pair(args):
    row, enzyme = args
    query = row['query']
    reference = row['reference']
    result = run_pair(query, reference, enzyme)
    rec = {
        'query': query,
        'reference': reference,
        'q_species': row.get('q_species'),
        'r_species': row.get('r_species'),
        'q_genus': row.get('q_genus'),
        'r_genus': row.get('r_genus'),
        'enzyme': enzyme,
    }
    if 'error' in result:
        rec['error'] = result['error']
    else:
        rec.update(result)
    return rec


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pairs', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--threads', type=int, default=min(16, cpu_count()))
    args = parser.parse_args()

    pairs = pd.read_csv(args.pairs, sep='\t')
    tasks = []
    for _, row in pairs.iterrows():
        for enzyme in ENZYMES:
            tasks.append((row.to_dict(), enzyme))

    total = len(tasks)
    print(f'Total tasks: {total}, threads: {args.threads}', flush=True)

    records = []
    with Pool(processes=args.threads) as pool:
        for i, rec in enumerate(pool.imap_unordered(process_pair, tasks), 1):
            records.append(rec)
            if i % 16 == 0 or i == total:
                print(f'[{i}/{total}] completed', flush=True)

    df = pd.DataFrame(records)
    df.to_csv(args.output, sep='\t', index=False, float_format='%.6f')
    print(f'Wrote {len(df)} records to {args.output}')


if __name__ == '__main__':
    main()
