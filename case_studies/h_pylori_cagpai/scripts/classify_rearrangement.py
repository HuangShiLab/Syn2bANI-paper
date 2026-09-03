#!/usr/bin/env python3
"""Extend cagPAI status with structural rearrangement calls from syn2bani struct --bed.

Optionally filters out genome-spanning calls that are likely circular-origin
artifacts (different arbitrary start coordinate on a circular chromosome).
"""
import argparse
import csv
import sys
from pathlib import Path

CAG_START = 547327
CAG_END = 583481
CAG_BUFFER = 2000  # include flanking SVs that may affect island boundaries
REGION = (CAG_START - CAG_BUFFER, CAG_END + CAG_BUFFER)


def overlaps_cagpai(chrom, start, end):
    return chrom == 'NC_000915.1' and start < REGION[1] and end > REGION[0]


def parse_bed(path, ref_length=None, artifact_threshold=None):
    """Parse BED; optionally exclude genome-spanning artifact calls.

    A call spanning > artifact_threshold fraction of ref_length is treated as a
    circular-origin coordinate artifact and excluded from the returned list (but
    counted separately for reporting).
    """
    svs = []
    excluded = []
    if not path.exists():
        return svs, excluded
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            if len(cols) < 4:
                continue
            chrom, start, end, name = cols[0], int(cols[1]), int(cols[2]), cols[3]
            span = end - start
            if (ref_length is not None and artifact_threshold is not None and
                    span > artifact_threshold * ref_length):
                sv_type = name.split('_')[0]
                excluded.append({
                    'chrom': chrom, 'start': start, 'end': end,
                    'name': name, 'type': sv_type,
                    'span': span,
                    'score': cols[4] if len(cols) > 4 else '.',
                })
                continue
            if overlaps_cagpai(chrom, start, end):
                sv_type = name.split('_')[0]
                svs.append({
                    'chrom': chrom, 'start': start, 'end': end,
                    'name': name, 'type': sv_type,
                    'score': cols[4] if len(cols) > 4 else '.',
                })
    return svs, excluded


def classify(status, svs):
    """Return extended state and summary string."""
    if status in ('partial', 'empty'):
        return status, svs
    # status == complete
    rearr = [s for s in svs if s['type'] in ('INV', 'TRA')]
    if rearr:
        return 'complete_rearranged', rearr
    large_del = [s for s in svs if s['type'] == 'DEL' and (s['end'] - s['start']) >= 10000]
    if large_del:
        # Presence/absence said complete, but struct sees a large deletion covering cagPAI.
        # Treat as partial (island integrity compromised).
        return 'partial', large_del
    return 'complete_collinear', []


def main():
    parser = argparse.ArgumentParser(
        description='Extend cagPAI marker status with syn2bani struct calls.'
    )
    parser.add_argument(
        '--states', type=Path,
        default=Path('/Volumes/MoneyCat/Data/song_2026_hpylori/cagpai_status/cagpai_states.tsv'),
        help='Input cagPAI marker status TSV.'
    )
    parser.add_argument(
        '--struct-dir', type=Path,
        default=Path('/Volumes/MoneyCat/Data/song_2026_hpylori/struct_vs_26695_filtered'),
        help='Directory containing <genome>.vs_hp26695.bed files.'
    )
    parser.add_argument(
        '--out', type=Path,
        default=Path('/Volumes/MoneyCat/Data/song_2026_hpylori/cagpai_status/cagpai_states_extended.tsv'),
        help='Output extended status TSV.'
    )
    parser.add_argument(
        '--filter-circular-origin', action='store_true',
        help='Exclude genome-spanning calls likely caused by circular-origin coordinate differences.'
    )
    parser.add_argument(
        '--ref-length', type=int, default=1667825,
        help='Reference chromosome length for artifact filtering (default: H. pylori 26695).'
    )
    parser.add_argument(
        '--artifact-threshold', type=float, default=0.50,
        help='Fraction of reference length above which a call is treated as artifact (default: 0.50).'
    )
    args = parser.parse_args()

    rows = []
    artifact_counts = {'excluded_calls': 0, 'flagged_genomes': 0}

    with open(args.states) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        fieldnames = reader.fieldnames + [
            'status_extended', 'n_cag_sv', 'cag_sv_types', 'cag_sv_list',
            'n_artifact_calls', 'artifact_call_list'
        ]
        for r in reader:
            gid = r['genome']
            bed = args.struct_dir / f'{gid}.vs_hp26695.bed'
            svs, excluded = parse_bed(
                bed,
                ref_length=args.ref_length if args.filter_circular_origin else None,
                artifact_threshold=args.artifact_threshold if args.filter_circular_origin else None,
            )
            artifact_counts['excluded_calls'] += len(excluded)
            if excluded:
                artifact_counts['flagged_genomes'] += 1

            ext, reported_svs = classify(r['status'], svs)
            r['status_extended'] = ext
            r['n_cag_sv'] = str(len(reported_svs))
            r['cag_sv_types'] = ';'.join(sorted(set(s['type'] for s in reported_svs))) if reported_svs else 'none'
            r['cag_sv_list'] = ';'.join(f"{s['type']}:{s['start']}-{s['end']}" for s in reported_svs) if reported_svs else 'none'
            r['n_artifact_calls'] = str(len(excluded))
            r['artifact_call_list'] = ';'.join(f"{s['type']}:{s['start']}-{s['end']}" for s in excluded) if excluded else 'none'
            rows.append(r)

    with open(args.out, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    # Print summary
    from collections import Counter
    cnt = Counter(r['status_extended'] for r in rows)
    print('Extended state counts:', file=sys.stderr)
    for st, n in sorted(cnt.items(), key=lambda x: -x[1]):
        print(f'  {st}: {n}', file=sys.stderr)
    if args.filter_circular_origin:
        print(f'Excluded {artifact_counts["excluded_calls"]} artifact calls '
              f'across {artifact_counts["flagged_genomes"]} genomes '
              f'(threshold {args.artifact_threshold:.0%} of {args.ref_length} bp).',
              file=sys.stderr)
    print(f'Wrote {args.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
