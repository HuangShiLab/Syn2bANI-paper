#!/usr/bin/env python3
"""Extract Syn2bANI raw-features for a set of genome pairs."""
import argparse
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd


def resolve_genome(accession: str, genomes_dir: Path) -> Path:
    direct = genomes_dir / f"{accession}.fna"
    if direct.exists():
        return direct
    matches = list(genomes_dir.glob(f"*{accession}*.fna"))
    if matches:
        return matches[0]
    raise FileNotFoundError(f"Genome not found for {accession}")


def run_syn2bani(args: dict) -> dict:
    q, r = args['query'], args['reference']
    try:
        cmd = [
            args['syn2bani'], 'dist',
            str(args['q_path']), str(args['r_path']),
            '--raw-features', '--min-af', '0.0',
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        for line in result.stdout.strip().split('\n'):
            if not line.strip() or line.startswith('#') or line.startswith('query_file'):
                continue
            parts = line.split('\t')
            if len(parts) >= 14:
                try:
                    return {
                        'query': q, 'reference': r,
                        's2b_raw_ani': float(parts[4]),
                        's2b_mash_ani': float(parts[5]),
                        's2b_chained_kmer_ani': float(parts[6]),
                        's2b_af_q': float(parts[7]),
                        's2b_af_r': float(parts[8]),
                        's2b_shared_tags': int(parts[9]),
                        's2b_containment': float(parts[10]),
                        's2b_ref_gc': float(parts[12]),
                        's2b_corrected_ani': float(parts[13]),
                        's2b_error': None,
                    }
                except (ValueError, IndexError):
                    continue
        return {'query': q, 'reference': r, 's2b_error': 'parse_failed'}
    except Exception as e:
        return {'query': q, 'reference': r, 's2b_error': str(e)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pairs', required=True)
    parser.add_argument('--genomes', required=True)
    parser.add_argument('--syn2bani', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--threads', type=int, default=8)
    args = parser.parse_args()

    genomes_dir = Path(args.genomes)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    pairs = pd.read_csv(args.pairs, sep='\t', usecols=['query', 'reference'])
    print(f'Loaded {len(pairs)} pairs', file=sys.stderr)

    records = []
    for _, row in pairs.iterrows():
        q_path = resolve_genome(row['query'], genomes_dir)
        r_path = resolve_genome(row['reference'], genomes_dir)
        records.append({
            'query': row['query'],
            'reference': row['reference'],
            'q_path': q_path,
            'r_path': r_path,
            'syn2bani': args.syn2bani,
        })

    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        results = list(pool.map(run_syn2bani, records))

    df = pd.DataFrame(results)
    df.to_csv(output_path, sep='\t', index=False, float_format='%.6f')
    print(f'Wrote {len(df)} rows to {output_path}', file=sys.stderr)
    errors = df['s2b_error'].notna().sum()
    if errors:
        print(f'Errors: {errors}', file=sys.stderr)


if __name__ == '__main__':
    main()
