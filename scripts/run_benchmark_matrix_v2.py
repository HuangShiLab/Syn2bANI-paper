#!/usr/bin/env python3
"""
run_benchmark_matrix_v2.py
Run skani, FastANI (pyfastani fallback), and Syn2bANI (raw + GBRT) on a set of genome pairs.

Supports checkpoint/resume: if --output exists, already-completed pairs are skipped.

Usage:
  python3 run_benchmark_matrix_v2.py \
    --pairs results/pairs_gtdb_r207.tsv \
    --genomes ~/data/gtdb-r207/genomes_all/ \
    --syn2bani ~/Downloads/Syn2bANI/target/release/syn2bani \
    --skani $(which skani) \
    --output results/matrix.tsv \
    --threads 16
"""

import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from multiprocessing import Pool
from pathlib import Path
from typing import Optional

import pandas as pd
import numpy as np
import shutil

try:
    import pyfastani
    HAS_PYFASTANI = True
except ImportError:
    HAS_PYFASTANI = False


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

def resolve_genome_path(accession: str, genomes_dir) -> Optional[Path]:
    """Resolve an accession to a genome file in genomes_dir."""
    genomes_dir = Path(genomes_dir)
    direct = genomes_dir / f"{accession}.fna"
    if direct.exists():
        return direct
    matches = list(genomes_dir.glob(f"*{accession}*.fna"))
    if matches:
        return matches[0]
    return None


# ---------------------------------------------------------------------------
# FastANI reference (pyfastani)
# ---------------------------------------------------------------------------


def _parse_fasta(path: Path) -> list[str]:
    """Parse a FASTA file into a list of contig sequences."""
    contigs = []
    current = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current:
                    contigs.append(''.join(current))
                current = []
            elif line:
                current.append(line)
        if current:
            contigs.append(''.join(current))
    return contigs


def run_fastani_pair_pyfastani(args: dict) -> dict:
    """Run FastANI reference via pyfastani."""
    q, r = args['query'], args['reference']
    q_path = resolve_genome_path(q, args['genomes_dir'])
    r_path = resolve_genome_path(r, args['genomes_dir'])

    if q_path is None or r_path is None:
        return {'query': q, 'reference': r, 'fastani_ani': None,
                'error': 'missing_genome'}

    try:
        # NOTE: do not cache pyfastani Sketch objects across queries;
        # reusing a sketch/index causes incorrect no-hit results.
        sketch = pyfastani.Sketch()
        sketch.add_draft('ref', _parse_fasta(r_path))
        mapper = sketch.index()
        hits = mapper.query_draft(_parse_fasta(q_path), threads=1)
        if hits:
            best = max(hits, key=lambda h: h.identity)
            return {'query': q, 'reference': r, 'fastani_ani': best.identity / 100.0,
                    'error': None}
        return {'query': q, 'reference': r, 'fastani_ani': 0.0,
                'error': 'no_alignment'}
    except Exception as e:
        return {'query': q, 'reference': r, 'fastani_ani': None,
                'error': str(e)}


def run_fastani_pair_binary(args: dict) -> dict:
    """Run FastANI binary if available."""
    q, r, fastani_path = args['query'], args['reference'], args['fastani_path']
    q_path = resolve_genome_path(q, args['genomes_dir'])
    r_path = resolve_genome_path(r, args['genomes_dir'])

    if q_path is None or r_path is None:
        return {'query': q, 'reference': r, 'fastani_ani': None,
                'error': 'missing_genome'}

    try:
        result = subprocess.run(
            [fastani_path, '-q', str(q_path), '-r', str(r_path),
             '-o', '/dev/stdout'],
            capture_output=True, text=True, timeout=300
        )
        lines = result.stdout.strip().split('\n')
        for line in lines:
            if not line.strip():
                continue
            parts = line.split('\t')
            if len(parts) >= 3:
                return {
                    'query': q, 'reference': r,
                    'fastani_ani': float(parts[2]) / 100.0,  # FastANI reports 0-100, normalize to 0-1
                    'error': None
                }
        return {'query': q, 'reference': r, 'fastani_ani': None,
                'error': 'parse_failed'}
    except Exception as e:
        return {'query': q, 'reference': r, 'fastani_ani': None,
                'error': str(e)}


# ---------------------------------------------------------------------------
# skani
# ---------------------------------------------------------------------------

def run_skani_pair(args: dict) -> dict:
    """Run skani dist on one pair."""
    q, r, skani_path = args['query'], args['reference'], args['skani_path']
    q_path = resolve_genome_path(q, args['genomes_dir'])
    r_path = resolve_genome_path(r, args['genomes_dir'])

    if q_path is None or r_path is None:
        return {'query': q, 'reference': r, 'skani_ani': None,
                'skani_align_frac': None, 'error': 'missing_genome'}

    try:
        result = subprocess.run(
            [skani_path, 'dist', str(q_path), str(r_path)],
            capture_output=True, text=True, timeout=120
        )
        lines = result.stdout.strip().split('\n')
        data_lines = []
        for line in lines:
            if not line.strip() or line.startswith('[') or 'INFO' in line:
                continue
            parts = line.split('\t')
            if len(parts) >= 3:
                try:
                    ani = float(parts[2])
                    align_frac = float(parts[3]) if len(parts) > 3 else None
                    return {
                        'query': q, 'reference': r,
                        'skani_ani': ani / 100.0,  # skani reports 0-100, normalize to 0-1
                        'skani_align_frac': align_frac,
                        'error': None
                    }
                except ValueError:
                    # Likely the header line; skip it
                    continue
        # If we only saw the header and no data line, skani found no alignment
        return {'query': q, 'reference': r, 'skani_ani': 0.0,
                'skani_align_frac': 0.0, 'error': None}
    except Exception as e:
        return {'query': q, 'reference': r, 'skani_ani': None,
                'skani_align_frac': None, 'error': str(e)}


# ---------------------------------------------------------------------------
# Syn2bANI
# ---------------------------------------------------------------------------

# TSV column indices from syn2bani dist output:
# 0 query_file, 1 ref_file, 2 ani, 3 af_q, 4 af_r,
# 5 query_name, 6 ref_name, 7 shared_tags, 8 sv_count
S2B_COL_ANI = 2
S2B_COL_AF_Q = 3
S2B_COL_AF_R = 4
S2B_COL_SHARED = 7


# raw-features TSV column indices:
# 0 query_file, 1 ref_file, 2 query_name, 3 ref_name,
# 4 raw_ani, 5 mash_ani, 6 chained_kmer_ani, 7 af_q, 8 af_r, 9 shared_tags,
# 10 containment, 11 div_proxy, 12 ref_gc, 13 corrected_ani
S2B_RAW_COL_ANI = 4
S2B_RAW_COL_MASH_ANI = 5
S2B_RAW_COL_CHAINED_KMER_ANI = 6
S2B_RAW_COL_AF_Q = 7
S2B_RAW_COL_AF_R = 8
S2B_RAW_COL_SHARED = 9
S2B_RAW_COL_REF_GC = 12
S2B_RAW_COL_CORRECTED = 13


def _parse_syn2bani_raw_features(stdout: str) -> dict:
    for line in stdout.strip().split('\n'):
        if line.startswith('#') or not line.strip():
            continue
        parts = line.split('\t')
        if parts[S2B_RAW_COL_ANI] == 'raw_ani':
            continue
        if len(parts) >= 14:
            return {
                'raw_ani': float(parts[S2B_RAW_COL_ANI]),
                'mash_ani': float(parts[S2B_RAW_COL_MASH_ANI]),
                'chained_kmer_ani': float(parts[S2B_RAW_COL_CHAINED_KMER_ANI]),
                'af_q': float(parts[S2B_RAW_COL_AF_Q]),
                'af_r': float(parts[S2B_RAW_COL_AF_R]),
                'shared_tags': int(parts[S2B_RAW_COL_SHARED]),
                'ref_gc': float(parts[S2B_RAW_COL_REF_GC]),
                'corrected_ani': float(parts[S2B_RAW_COL_CORRECTED]),
            }
    return None


def run_syn2bani_pair(args: dict) -> dict:
    """Run Syn2bANI dist --raw-features on one pair (raw + GBRT in one shot)."""
    q, r, syn2bani_path = args['query'], args['reference'], args['syn2bani_path']
    q_path = resolve_genome_path(q, args['genomes_dir'])
    r_path = resolve_genome_path(r, args['genomes_dir'])

    if q_path is None or r_path is None:
        return {'query': q, 'reference': r, 's2b_raw_ani': None,
                's2b_gbrt_ani': None, 's2b_shared_tags': None,
                'error': 'missing_genome'}

    results = {'query': q, 'reference': r}

    try:
        cmd = [syn2bani_path, 'dist', str(q_path), str(r_path),
               '-e', 'AloI,BslFI', '--raw-features', '--min-af', '0.0']
        result = subprocess.run(
            cmd,
            capture_output=True, text=True, timeout=120
        )
        parsed = _parse_syn2bani_raw_features(result.stdout)
        if parsed:
            results['s2b_raw_ani'] = parsed['raw_ani']
            results['s2b_mash_ani'] = parsed['mash_ani']
            results['s2b_chained_kmer_ani'] = parsed['chained_kmer_ani']
            results['s2b_gbrt_ani'] = parsed['corrected_ani']
            results['s2b_shared_tags'] = parsed['shared_tags']
            results['s2b_af_q'] = parsed['af_q']
            results['s2b_af_r'] = parsed['af_r']
            results['s2b_ref_gc'] = parsed['ref_gc']
        else:
            results['s2b_raw_ani'] = None
            results['s2b_gbrt_ani'] = None
            results['s2b_shared_tags'] = None
            results['error'] = f'parse_failed rc={result.returncode} stdout={result.stdout[:200]} stderr={result.stderr[:200]}'
    except Exception as e:
        results['s2b_raw_ani'] = None
        results['s2b_gbrt_ani'] = None
        results['s2b_shared_tags'] = None
        results['error'] = str(e)

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pairs', required=True)
    parser.add_argument('--genomes', required=True)
    parser.add_argument('--syn2bani', default='syn2bani')
    parser.add_argument('--skani', default='skani')
    parser.add_argument('--fastani', default='fastANI')
    parser.add_argument('--use-pyfastani', action='store_true',
                        help='Use pyfastani Python package instead of fastANI binary')
    parser.add_argument('--output', required=True)
    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--tools', default='all',
                        help='Comma-separated: skani,fastani,syn2bani,all')
    parser.add_argument('--chunk-size', type=int, default=100,
                        help='Save checkpoint every N pairs')
    args = parser.parse_args()

    genomes_dir = Path(args.genomes)
    if not genomes_dir.exists():
        print(f"ERROR: Genomes directory not found: {genomes_dir}")
        sys.exit(1)

    pairs_df = pd.read_csv(args.pairs, sep='\t')
    print(f"Loaded {len(pairs_df)} pairs from {args.pairs}")

    # Resume from existing output
    output_path = Path(args.output)
    completed = set()
    existing = None
    if output_path.exists():
        existing = pd.read_csv(output_path, sep='\t')
        for _, row in existing.iterrows():
            completed.add((row['query'], row['reference']))
        print(f"Resuming: {len(completed)} pairs already completed")

    remaining = pairs_df[~pairs_df.apply(
        lambda row: (row['query'], row['reference']) in completed, axis=1
    )].copy()
    print(f"Remaining pairs to run: {len(remaining)}")

    if remaining.empty:
        print("All pairs already completed.")
        return

    tools = [t.strip() for t in args.tools.split(',')]
    run_all = 'all' in tools

    # --- Run FastANI ---
    if run_all or 'fastani' in tools:
        use_pyfastani = args.use_pyfastani or not shutil.which(args.fastani)
        if use_pyfastani and not HAS_PYFASTANI:
            print("WARNING: pyfastani not available and fastANI binary not found; skipping FastANI")
        else:
            print(f"\n[1/3] Running FastANI on {len(remaining)} pairs...")
            runner = run_fastani_pair_pyfastani if use_pyfastani else run_fastani_pair_binary
            fastani_args = [{
                'query': row['query'], 'reference': row['reference'],
                'genomes_dir': genomes_dir, 'fastani_path': args.fastani
            } for _, row in remaining.iterrows()]

            # Use parallel pool for FastANI binary (process-safe on Linux);
            # keep sequential fallback for pyfastani on macOS.
            if use_pyfastani:
                print("  (pyfastani runs sequentially to avoid concurrency issues)")
                fastani_results = [runner(a) for a in fastani_args]
            else:
                print(f"  (Running FastANI binary in parallel with {args.threads} processes)")
                with Pool(args.threads) as pool:
                    fastani_results = pool.map(runner, fastani_args)

            fastani_df = pd.DataFrame(fastani_results)
            remaining = remaining.merge(
                fastani_df[['query', 'reference', 'fastani_ani']],
                on=['query', 'reference'], how='left'
            )
            print(f"  FastANI completed. Missing: {remaining['fastani_ani'].isna().sum()}")

    # --- Run skani ---
    if run_all or 'skani' in tools:
        print(f"\n[2/3] Running skani on {len(remaining)} pairs...")
        skani_args = [{
            'query': row['query'], 'reference': row['reference'],
            'genomes_dir': genomes_dir, 'skani_path': args.skani
        } for _, row in remaining.iterrows()]

        with Pool(args.threads) as pool:
            skani_results = pool.map(run_skani_pair, skani_args)

        skani_df = pd.DataFrame(skani_results)
        remaining = remaining.merge(
            skani_df[['query', 'reference', 'skani_ani', 'skani_align_frac']],
            on=['query', 'reference'], how='left'
        )
        print(f"  skani completed. Missing: {remaining['skani_ani'].isna().sum()}")

    # --- Run Syn2bANI ---
    if run_all or 'syn2bani' in tools:
        print(f"\n[3/3] Running Syn2bANI on {len(remaining)} pairs...")
        s2b_args = [{
            'query': row['query'], 'reference': row['reference'],
            'genomes_dir': genomes_dir, 'syn2bani_path': args.syn2bani
        } for _, row in remaining.iterrows()]

        with Pool(args.threads) as pool:
            s2b_results = pool.map(run_syn2bani_pair, s2b_args)

        s2b_df = pd.DataFrame(s2b_results)
        remaining = remaining.merge(
            s2b_df[['query', 'reference', 's2b_raw_ani', 's2b_mash_ani', 's2b_chained_kmer_ani',
                    's2b_gbrt_ani', 's2b_shared_tags', 's2b_af_q', 's2b_af_r', 's2b_ref_gc']],
            on=['query', 'reference'], how='left'
        )
        print(f"  Syn2bANI completed. Missing raw: {remaining['s2b_raw_ani'].isna().sum()}")

    # --- Merge with existing and save ---
    os.makedirs(output_path.parent or '.', exist_ok=True)
    if output_path.exists():
        combined = pd.concat([existing, remaining], ignore_index=True)
    else:
        combined = remaining

    combined.to_csv(args.output, sep='\t', index=False, float_format='%.6f')
    print(f"\n{'='*60}")
    print(f"Results saved: {args.output}")
    print(f"Total rows: {len(combined)}")
    print(f"Columns: {list(combined.columns)}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
