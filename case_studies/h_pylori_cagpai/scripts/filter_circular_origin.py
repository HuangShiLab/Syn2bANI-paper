#!/usr/bin/env python3
"""Flag H. pylori structural calls that are likely circular-origin artifacts.

A call spanning more than a fixed fraction of the reference chromosome is
interpreted as a coordinate-system difference (different arbitrary start on a
circular chromosome) rather than a biological rearrangement. Such calls overlap
the cagPAI window by construction and inflate the `complete_rearranged` state.

Input: BED files produced by `syn2bani struct --bed` (one per query genome).
Output: TSV with columns genome, max_span, ref_length, max_span_fraction,
        flagged_artifact, largest_call_type, largest_call_coords.
"""

import argparse
import csv
import sys
from pathlib import Path


def parse_bed(path):
    """Return list of (chrom, start, end, name) tuples from a BED file."""
    svs = []
    if not path.exists():
        return svs
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 4:
                continue
            chrom, start, end, name = cols[0], int(cols[1]), int(cols[2]), cols[3]
            svs.append((chrom, start, end, name))
    return svs


def main():
    parser = argparse.ArgumentParser(
        description='Flag circular-origin artifacts in syn2bani struct BED output.'
    )
    parser.add_argument(
        '--struct-dir', required=True, type=Path,
        help='Directory containing <genome>.vs_hp26695.bed files.'
    )
    parser.add_argument(
        '--ref-length', type=int, default=1667825,
        help='Length of the reference chromosome (default: H. pylori 26695, 1,667,825 bp).'
    )
    parser.add_argument(
        '--threshold', type=float, default=0.50,
        help='Fraction of reference length above which a call is flagged (default: 0.50).'
    )
    parser.add_argument(
        '--out', type=Path, default=Path('circular_origin_flags.tsv'),
        help='Output TSV path.'
    )
    args = parser.parse_args()

    if not args.struct_dir.is_dir():
        print(f'Error: --struct-dir does not exist: {args.struct_dir}', file=sys.stderr)
        sys.exit(1)

    fieldnames = [
        'genome', 'n_sv', 'max_span', 'ref_length', 'max_span_fraction',
        'flagged_artifact', 'largest_call_type', 'largest_call_coords'
    ]
    rows = []

    for bed_path in sorted(args.struct_dir.glob('*.vs_hp26695.bed')):
        genome = bed_path.stem.replace('.vs_hp26695', '')
        svs = parse_bed(bed_path)
        if not svs:
            rows.append({
                'genome': genome,
                'n_sv': 0,
                'max_span': 0,
                'ref_length': args.ref_length,
                'max_span_fraction': 0.0,
                'flagged_artifact': False,
                'largest_call_type': 'none',
                'largest_call_coords': 'none',
            })
            continue

        spans = []
        for chrom, start, end, name in svs:
            span = end - start
            sv_type = name.split('_')[0]
            spans.append((span, sv_type, f'{chrom}:{start}-{end}'))

        spans.sort(reverse=True)
        max_span, max_type, max_coords = spans[0]
        fraction = max_span / args.ref_length
        flagged = fraction > args.threshold

        rows.append({
            'genome': genome,
            'n_sv': len(svs),
            'max_span': max_span,
            'ref_length': args.ref_length,
            'max_span_fraction': round(fraction, 4),
            'flagged_artifact': flagged,
            'largest_call_type': max_type,
            'largest_call_coords': max_coords,
        })

    with open(args.out, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    n_flagged = sum(1 for r in rows if r['flagged_artifact'])
    print(f'Wrote {args.out}', file=sys.stderr)
    print(f'Genomes: {len(rows)}; flagged (max span > {args.threshold:.0%} of ref): {n_flagged}',
          file=sys.stderr)


if __name__ == '__main__':
    main()
