#!/usr/bin/env python3
"""Uniform-sim gate: the candidate models must not degrade the simindel ladder
(MAE must stay ~0.07). Flexible mixtures can invent heterogeneity on uniform
data, so this gate decides whether a parsimony gate (LRT vs gamma) is needed.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from explore_mosaic import load_strata
from model_mosaic import kernel_tensor, counts_from_strata, npmle_em, \
    load_ani_row, A_GRID
import model_mosaic as mm
from model_discrete import discrete_fit

HERE = os.path.dirname(os.path.abspath(__file__))
SD = os.path.join(HERE, "simindel_strata")
MANIFEST = "/Users/macstudio/Downloads/Syn2bANI/prototype/simindel/manifest.tsv"


def main():
    truth = {}
    with open(MANIFEST) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            truth[p[0]] = float(p[1]) * 100

    # cap the grid once for all NPcap fits
    mm.V_GRID = mm.V_GRID[mm.V_GRID <= 0.450001]
    mm.A_GRID = np.exp(-mm.V_GRID)

    print(f"{'case':<16} {'truth':>6} {'gamma':>7} {'unif':>7} {'D2':>7} {'NPcap':>7}")
    sums = {}
    for name in sorted(truth):
        strata = load_strata(os.path.join(SD, name + ".strata.tsv"))
        row = load_ani_row(os.path.join(SD, name + ".ani.tsv"))
        counts = counts_from_strata(strata)
        K = kernel_tensor(strata)
        vs, ws, _ = discrete_fit(counts, K, 2)
        d2 = 100 * float(ws @ np.exp(-vs))
        wnp = npmle_em(counts, K, iters=5000)
        npc = 100 * float(wnp @ mm.A_GRID)
        g, u = float(row["ani"]), float(row["ani_uniform"])
        print(f"{name:<16} {truth[name]:6.2f} {g:7.3f} {u:7.3f} {d2:7.3f} {npc:7.3f}")
        for k, v in [("gamma", g), ("uniform", u), ("D2", d2), ("NPcap", npc)]:
            sums.setdefault(k, []).append(v - truth[name])
    print()
    for k, errs in sums.items():
        errs = np.array(errs)
        print(f"{k:<8} MAE {np.abs(errs).mean():6.4f}  bias {errs.mean():+6.4f}")


if __name__ == "__main__":
    main()
