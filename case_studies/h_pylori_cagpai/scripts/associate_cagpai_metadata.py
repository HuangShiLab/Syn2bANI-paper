#!/usr/bin/env python3
"""Associate extended cagPAI states with H. pylori isolate metadata."""
import csv
import sys
from pathlib import Path
from collections import Counter
import numpy as np

BASE = Path('/Volumes/MoneyCat/Data/song_2026_hpylori/cagpai_status')
META = Path('/Volumes/MoneyCat/Data/song_2026_hpylori/metadata.csv')
STATES = BASE / 'cagpai_states_extended.tsv'
OUT = BASE / 'cagpai_association.tsv'
FIGDIR = BASE / 'figures'
FIGDIR.mkdir(parents=True, exist_ok=True)

STATE_ORDER = ['empty', 'partial', 'complete_collinear', 'complete_rearranged']
STATE_COLORS = {
    'empty': '#d62728',
    'partial': '#ff7f0e',
    'complete_collinear': '#2ca02c',
    'complete_rearranged': '#1f77b4',
}


def load_states():
    d = {}
    with open(STATES) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            d[r['genome']] = r
    return d


def load_metadata():
    rows = []
    with open(META) as fh:
        for r in csv.DictReader(fh):
            rows.append(r)
    return rows


def chisq(table):
    """table: dict group->Counter(state). Returns chi2, p, dof."""
    try:
        from scipy.stats import chi2_contingency
    except ImportError:
        return None, None, None
    groups = sorted(table.keys())
    arr = np.array([[table[g].get(s, 0) for s in STATE_ORDER] for g in groups])
    if arr.shape[0] < 2 or arr.shape[1] < 2 or (arr.sum(axis=1) < 5).any():
        return None, None, None
    chi2, p, dof, _ = chi2_contingency(arr)
    return chi2, p, dof


def summarize(states, meta, field):
    table = {}
    missing = 0
    for r in meta:
        gid = r['assembly']
        val = r.get(field, '')
        if not val:
            missing += 1
            continue
        st = states.get(gid, {}).get('status_extended', 'unknown')
        table.setdefault(val, Counter())[st] += 1
    return table, missing


def write_table(fh, field, table):
    fh.write(f'\n## {field}\n')
    groups = sorted(table.keys())
    fh.write('group\t' + '\t'.join(STATE_ORDER) + '\ttotal\n')
    for g in groups:
        cnt = table[g]
        total = sum(cnt.values())
        frac = [cnt.get(s, 0) / total for s in STATE_ORDER]
        fh.write(f'{g}\t' + '\t'.join(f'{cnt.get(s,0)} ({f:.2%})' for s, f in zip(STATE_ORDER, frac)) + f'\t{total}\n')
    chi2, p, dof = chisq(table)
    if chi2 is not None:
        fh.write(f'chi2={chi2:.3f}, df={dof}, p={p:.4g}\n')
    else:
        fh.write('chi2=NA (insufficient counts or scipy missing)\n')


def plot_stacked(table, field, outpath):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib not available, skipping plot', file=sys.stderr)
        return
    groups = sorted(table.keys())
    arr = np.array([[table[g].get(s, 0) for s in STATE_ORDER] for g in groups], dtype=float)
    totals = arr.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1
    frac = arr / totals

    fig, ax = plt.subplots(figsize=(max(6, len(groups) * 0.6), 5))
    left = np.zeros(len(groups))
    for i, st in enumerate(STATE_ORDER):
        ax.barh(groups, frac[:, i], left=left, color=STATE_COLORS[st], label=st)
        left += frac[:, i]
    ax.set_xlim(0, 1)
    ax.set_xlabel('proportion')
    ax.set_title(f'cagPAI extended state by {field}')
    ax.legend(loc='lower right', title='state')
    plt.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def main():
    states = load_states()
    meta = load_metadata()

    with open(OUT, 'w') as fh:
        fh.write(f'# cagPAI state–metadata associations (n={len(meta)})\n')
        overall = Counter(states.get(r['assembly'], {}).get('status_extended', 'unknown') for r in meta)
        fh.write('\n## overall\n')
        for s in STATE_ORDER:
            fh.write(f'{s}\t{overall.get(s, 0)} ({overall.get(s, 0) / len(meta):.2%})\n')

        for field in ['group', 'country', 'fastbaps', 'phylogenetic_population']:
            table, missing = summarize(states, meta, field)
            if len(table) < 2:
                continue
            write_table(fh, field, table)
            plot_stacked(table, field, FIGDIR / f'cagpai_state_by_{field}.png')

    print(f'Wrote {OUT}', file=sys.stderr)
    print(f'Figures in {FIGDIR}', file=sys.stderr)


if __name__ == '__main__':
    main()
