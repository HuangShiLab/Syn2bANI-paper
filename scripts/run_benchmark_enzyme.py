#!/usr/bin/env python3
"""Run Syn2bANI with a specific enzyme on a set of pairs and output raw-features."""
import argparse
import subprocess
from pathlib import Path
from multiprocessing import Pool

import pandas as pd


def resolve_genome_path(accession, genomes_dir):
    genomes_dir = Path(genomes_dir)
    direct = genomes_dir / f"{accession}.fna"
    if direct.exists():
        return direct
    matches = list(genomes_dir.glob(f"*{accession}*.fna"))
    if matches:
        return matches[0]
    return None


def run_syn2bani_pair(args):
    q, r, genomes_dir, syn2bani_path, enzyme = args['query'], args['reference'], args['genomes_dir'], args['syn2bani_path'], args['enzyme']
    q_path = resolve_genome_path(q, genomes_dir)
    r_path = resolve_genome_path(r, genomes_dir)

    if q_path is None or r_path is None:
        return {'query': q, 'reference': r, 'error': 'missing_genome'}

    cmd = [syn2bani_path, 'dist', str(q_path), str(r_path),
           '-e', enzyme, '--raw-features', '--min-af', '0.0']
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            return {'query': q, 'reference': r, 'error': result.stderr[:200]}

        # Parse raw-features TSV
        # 0 query_file, 1 ref_file, 2 query_name, 3 ref_name,
        # 4 raw_ani, 5 mash_ani, 6 af_q, 7 af_r, 8 shared_tags, 9 containment,
        # 10 div_proxy, 11 ref_gc, 12 corrected_ani
        for line in result.stdout.strip().split('\n'):
            if line.startswith('#') or not line.strip() or line.startswith('query_file'):
                continue
            parts = line.split('\t')
            if len(parts) >= 13:
                return {
                    'query': q, 'reference': r,
                    'raw_ani': float(parts[4]),
                    'mash_ani': float(parts[5]),
                    'af_q': float(parts[6]),
                    'af_r': float(parts[7]),
                    'shared_tags': int(parts[8]),
                    'ref_gc': float(parts[11]),
                    'corrected_ani': float(parts[12]),
                }
        return {'query': q, 'reference': r, 'error': 'parse_failed'}
    except Exception as e:
        return {'query': q, 'reference': r, 'error': str(e)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pairs', required=True)
    parser.add_argument('--genomes', required=True)
    parser.add_argument('--syn2bani', default='/lustre1/g/aos_shihuang/Syn2bANI/target/release/syn2bani')
    parser.add_argument('--enzyme', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--threads', type=int, default=16)
    args = parser.parse_args()

    pairs = pd.read_csv(args.pairs, sep='\t')
    print(f"Loaded {len(pairs)} pairs")

    tasks = [{
        'query': row['query'], 'reference': row['reference'],
        'genomes_dir': args.genomes, 'syn2bani_path': args.syn2bani,
        'enzyme': args.enzyme
    } for _, row in pairs.iterrows()]

    with Pool(args.threads) as pool:
        results = pool.map(run_syn2bani_pair, tasks)

    df = pd.DataFrame(results)
    # Prefix columns with enzyme name to avoid conflicts when merging
    prefix = args.enzyme.lower()
    rename = {}
    for col in ['raw_ani', 'mash_ani', 'af_q', 'af_r', 'shared_tags', 'ref_gc', 'corrected_ani']:
        if col in df.columns:
            rename[col] = f'{prefix}_{col}'
    df = df.rename(columns=rename)

    df.to_csv(args.output, sep='\t', index=False, float_format='%.6f')
    print(f"Wrote {len(df)} records to {args.output}")


if __name__ == '__main__':
    main()
