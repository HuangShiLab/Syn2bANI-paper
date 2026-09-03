#!/usr/bin/env python3
"""Lineage-stratified association between cagPAI state and disease stage.

The raw disease-stage association is confounded by FastBAPS lineage, which is
itself associated with both geography and cagPAI architecture. This script
computes Cochran-Mantel-Haenszel (CMH) statistics for 2x2 contrasts stratified
by FastBAPS lineage, where the conditioning removes the lineage-driven component.
"""
import csv
import sys
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np

try:
    from scipy.stats import chi2
except ImportError:
    chi2 = None

BASE = Path('/Volumes/MoneyCat/Data/song_2026_hpylori/cagpai_status')
META = Path('/Volumes/MoneyCat/Data/song_2026_hpylori/metadata.csv')
STATES = BASE / 'cagpai_states_extended.tsv'
OUT = BASE / 'cagpai_association_stratified.tsv'


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


def cmh_2x2xk(tables):
    """Compute Cochran-Mantel-Haenszel for a list of 2x2 tables.

    Each table is (a, b, c, d):
        a = case + outcome,    b = case + no outcome
        c = control + outcome, d = control + no outcome

    Returns chi2, p, dof=1, common odds ratio.
    """
    if chi2 is None:
        return None, None, 1, None

    sum_a = 0.0
    sum_ea = 0.0
    sum_var = 0.0
    or_num = 0.0
    or_den = 0.0

    for a, b, c, d in tables:
        n = a + b + c + d
        if n == 0:
            continue
        n1 = a + b      # case total
        n0 = c + d      # control total
        m1 = a + c      # outcome total
        m0 = b + d      # no-outcome total

        ea = n1 * m1 / n
        va = n1 * n0 * m1 * m0 / (n * n * (n - 1)) if n > 1 else 0.0

        sum_a += a
        sum_ea += ea
        sum_var += va

        # Mantel-Haenszel OR contribution
        or_den += b * c / n
        or_num += a * d / n

    if sum_var == 0:
        return None, None, 1, None

    chi2_stat = (abs(sum_a - sum_ea) - 0.5) ** 2 / sum_var  # continuity correction
    p = chi2.sf(chi2_stat, 1)
    or_mh = or_num / or_den if or_den > 0 else float('nan')

    return chi2_stat, p, 1, or_mh


def build_tables(rows, states, case_stages, control_stages, outcome_fn):
    """Build 2x2 tables per FastBAPS lineage."""
    strata = defaultdict(lambda: [[0, 0], [0, 0]])  # [case/control][outcome/no]
    for r in rows:
        gid = r['assembly']
        stage = r.get('group', '')
        lineage = r.get('fastbaps', '')
        if not stage or not lineage:
            continue
        if stage in case_stages:
            case_idx = 0
        elif stage in control_stages:
            case_idx = 1
        else:
            continue
        st = states.get(gid, {}).get('status_extended', 'unknown')
        if st == 'unknown':
            continue
        outcome = outcome_fn(st)
        outcome_idx = 0 if outcome else 1
        strata[lineage][case_idx][outcome_idx] += 1
    return strata


def print_contrast(fh, name, rows, states, case_stages, control_stages, outcome_fn):
    strata = build_tables(rows, states, case_stages, control_stages, outcome_fn)
    tables = []
    fh.write(f'\n### {name}: {"/".join(case_stages)} vs {"/".join(control_stages)}\n')
    fh.write('fastbaps\tcase_outcome\tcase_no_outcome\tcontrol_outcome\tcontrol_no_outcome\n')
    for lineage in sorted(strata):
        t = strata[lineage]
        a, b = t[0]  # case: outcome, no outcome
        c, d = t[1]  # control: outcome, no outcome
        fh.write(f'{lineage}\t{a}\t{b}\t{c}\t{d}\n')
        if all(v >= 1 for v in (a, b, c, d)):
            tables.append((a, b, c, d))
    if tables:
        chi2_stat, p, dof, or_mh = cmh_2x2xk(tables)
        fh.write(f'CMH chi2={chi2_stat:.4f}, df={dof}, p={p:.4g}, OR_MH={or_mh:.3f}\n')
    else:
        fh.write('CMH=NA (no stratum with all cells >= 1)\n')


def main():
    states = load_states()
    meta = load_metadata()

    with open(OUT, 'w') as fh:
        fh.write(f'# Lineage-stratified cagPAI–disease associations (n={len(meta)})\n')
        fh.write('# Extended state counts after circular-origin filtering:\n')
        overall = Counter(states.get(r['assembly'], {}).get('status_extended', 'unknown') for r in meta)
        for s in ['empty', 'partial', 'complete_collinear', 'complete_rearranged']:
            fh.write(f'#   {s}: {overall.get(s, 0)} ({overall.get(s, 0) / len(meta):.2%})\n')

        # Contrast 1: cagPAI presence (empty/partial = absent, complete_* = present)
        print_contrast(
            fh, 'cagPAI presence', meta, states, {'GC'}, {'NAG'},
            lambda st: st in ('complete_collinear', 'complete_rearranged')
        )

        # Contrast 2: cagPAI rearrangement among present (complete_rearranged vs complete_collinear)
        print_contrast(
            fh, 'cagPAI rearrangement', meta, states, {'GC'}, {'NAG'},
            lambda st: st == 'complete_rearranged'
        )

        # Contrast 3: advanced vs early for presence
        print_contrast(
            fh, 'cagPAI presence (advanced vs early)', meta, states, {'GC', 'IM'}, {'NAG', 'AG'},
            lambda st: st in ('complete_collinear', 'complete_rearranged')
        )

        # Contrast 4: advanced vs early for rearrangement
        print_contrast(
            fh, 'cagPAI rearrangement (advanced vs early)', meta, states, {'GC', 'IM'}, {'NAG', 'AG'},
            lambda st: st == 'complete_rearranged'
        )

    print(f'Wrote {OUT}', file=sys.stderr)


if __name__ == '__main__':
    main()
