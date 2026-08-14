#!/usr/bin/env python3
"""LRT gate between the gamma fit and the capped NPMLE, over all 31 sims
(19 mosaic + 12 uniform). The gate must separate the regimes: fire on mosaic
(real heterogeneity misspecification), not fire on uniform.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from explore_mosaic import load_strata, MOSAIC
from model_mosaic import kernel_tensor, counts_from_strata, npmle_em, \
    load_ani_row, collect_cases
import model_mosaic as mm
from het_fit import estimate_het, Stratum

HERE = os.path.dirname(os.path.abspath(__file__))
SD = os.path.join(HERE, "simindel_strata")
MANIFEST = "/Users/macstudio/Downloads/Syn2bANI/prototype/simindel/manifest.tsv"

mm.V_GRID = mm.V_GRID[mm.V_GRID <= 0.450001]
mm.A_GRID = np.exp(-mm.V_GRID)


def gamma_nll(counts, strata):
    """NLL of the data at the gamma MLE (supported or not)."""
    st = [Stratum(e, s["tag_len"], s["body_len"], list(map(int, s["hist"])),
                  int(s["n_miss"])) for e, s in strata.items()]
    h = estimate_het(st)
    return h["nll_het"], h


def npmle_nll(counts, K, w):
    ll = 0.0
    for e, c in counts.items():
        P = K[e] @ w
        ll += float(c @ np.log(np.maximum(P, 1e-300)))
    return -ll


def main():
    rows = []
    for c in collect_cases():
        strata = c["strata"]
        counts = counts_from_strata(strata)
        K = kernel_tensor(strata)
        w = npmle_em(counts, K, iters=5000)
        npc = 100 * float(w @ mm.A_GRID)
        nll_np = npmle_nll(counts, K, w)
        nll_g, h = gamma_nll(counts, strata)
        lrt = 2 * (nll_g - nll_np)
        rows.append(("mosaic", c["name"], c["truth"] * 100,
                     float(c["row"]["ani"]), npc, lrt, h["alpha_raw"]))

    truth = {}
    with open(MANIFEST) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            truth[p[0]] = float(p[1]) * 100
    for name in sorted(truth):
        strata = load_strata(os.path.join(SD, name + ".strata.tsv"))
        row = load_ani_row(os.path.join(SD, name + ".ani.tsv"))
        counts = counts_from_strata(strata)
        K = kernel_tensor(strata)
        w = npmle_em(counts, K, iters=5000)
        npc = 100 * float(w @ mm.A_GRID)
        nll_np = npmle_nll(counts, K, w)
        nll_g, h = gamma_nll(counts, strata)
        lrt = 2 * (nll_g - nll_np)
        rows.append(("uniform", name, truth[name], float(row["ani"]), npc, lrt,
                     h["alpha_raw"]))

    print(f"{'set':<8} {'case':<34} {'truth':>6} {'gamma':>7} {'NPcap':>7} {'LRT':>8}")
    for r in rows:
        print(f"{r[0]:<8} {r[1]:<34} {r[2]:6.2f} {r[3]:7.2f} {r[4]:7.2f} {r[5]:8.1f}")
    for thr in (5, 10, 20, 40, 80):
        for s in ("mosaic", "uniform"):
            sel = [r for r in rows if r[0] == s]
            est = [r[4] if r[5] > thr else r[3] for r in sel]
            errs = np.array([e - r[2] for e, r in zip(est, sel)])
            print(f"thr {thr:>3} {s:<8} MAE {np.abs(errs).mean():6.3f} "
                  f"bias {errs.mean():+6.3f}  fired {sum(1 for r in sel if r[5] > thr)}/{len(sel)}")


if __name__ == "__main__":
    main()
