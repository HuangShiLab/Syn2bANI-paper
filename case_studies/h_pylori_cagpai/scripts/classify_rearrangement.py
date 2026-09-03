#!/usr/bin/env python3
"""Extend cagPAI status with structural rearrangement calls from syn2bani struct --bed."""
import csv
import sys
from pathlib import Path

CAG_START = 547327
CAG_END = 583481
CAG_BUFFER = 2000  # include flanking SVs that may affect island boundaries
REGION = (CAG_START - CAG_BUFFER, CAG_END + CAG_BUFFER)

STATES = Path('/Volumes/MoneyCat/Data/song_2026_hpylori/cagpai_status/cagpai_states.tsv')
STRUCT = Path('/Volumes/MoneyCat/Data/song_2026_hpylori/struct_vs_26695')
OUT = Path('/Volumes/MoneyCat/Data/song_2026_hpylori/cagpai_status/cagpai_states_extended.tsv')


def overlaps_cagpai(chrom, start, end):
    return chrom == 'NC_000915.1' and start < REGION[1] and end > REGION[0]


def parse_bed(path):
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
            if overlaps_cagpai(chrom, start, end):
                sv_type = name.split('_')[0]
                svs.append({
                    'chrom': chrom, 'start': start, 'end': end,
                    'name': name, 'type': sv_type,
                    'score': cols[4] if len(cols) > 4 else '.',
                })
    return svs


def classify(status, svs):
    """Return extended state and summary string."""
    cag_svs = [s for s in svs]
    if status in ('partial', 'empty'):
        return status, cag_svs
    # status == complete
    rearr = [s for s in cag_svs if s['type'] in ('INV', 'TRA')]
    if rearr:
        return 'complete_rearranged', cag_svs
    large_del = [s for s in cag_svs if s['type'] == 'DEL' and (s['end'] - s['start']) >= 10000]
    if large_del:
        # Presence/absence said complete, but struct sees a large deletion covering cagPAI.
        # Treat as partial (island integrity compromised).
        return 'partial', cag_svs
    return 'complete_collinear', cag_svs


def main():
    rows = []
    with open(STATES) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        fieldnames = reader.fieldnames + ['status_extended', 'n_cag_sv', 'cag_sv_types', 'cag_sv_list']
        for r in reader:
            gid = r['genome']
            bed = STRUCT / f'{gid}.vs_hp26695.bed'
            svs = parse_bed(bed)
            ext, _ = classify(r['status'], svs)
            r['status_extended'] = ext
            r['n_cag_sv'] = str(len(svs))
            r['cag_sv_types'] = ';'.join(sorted(set(s['type'] for s in svs))) if svs else 'none'
            r['cag_sv_list'] = ';'.join(f"{s['type']}:{s['start']}-{s['end']}" for s in svs) if svs else 'none'
            rows.append(r)

    with open(OUT, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    # Print summary
    from collections import Counter
    cnt = Counter(r['status_extended'] for r in rows)
    print('Extended state counts:', file=sys.stderr)
    for st, n in sorted(cnt.items(), key=lambda x: -x[1]):
        print(f'  {st}: {n}', file=sys.stderr)
    print(f'Wrote {OUT}', file=sys.stderr)


if __name__ == '__main__':
    main()
