#!/usr/bin/env python3
"""
Run one chunk of the benchmark matrix on HPC.
Thin wrapper around run_benchmark_matrix_v2.py for job-array use.
"""
import argparse
import sys
from pathlib import Path

# Re-use the matrix runner logic
sys.path.insert(0, str(Path(__file__).parent))
from run_benchmark_matrix_v2 import main as matrix_main


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pairs', required=True, help='Chunk pair TSV')
    parser.add_argument('--genomes', required=True)
    parser.add_argument('--syn2bani', default='syn2bani')
    parser.add_argument('--skani', default='skani')
    parser.add_argument('--fastani', default='fastANI')
    parser.add_argument('--use-pyfastani', action='store_true')
    parser.add_argument('--output', required=True)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--tools', default='all')
    parser.add_argument('--chunk-size', type=int, default=1000)
    args = parser.parse_args()

    # Forward to the existing runner
    sys.argv = [
        'run_benchmark_matrix_v2.py',
        '--pairs', args.pairs,
        '--genomes', args.genomes,
        '--syn2bani', args.syn2bani,
        '--skani', args.skani,
        '--fastani', args.fastani,
        '--output', args.output,
        '--threads', str(args.threads),
        '--tools', args.tools,
        '--chunk-size', str(args.chunk_size),
    ]
    if args.use_pyfastani:
        sys.argv.append('--use-pyfastani')
    matrix_main()


if __name__ == '__main__':
    main()
