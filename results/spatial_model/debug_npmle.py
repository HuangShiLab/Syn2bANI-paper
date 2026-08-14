#!/usr/bin/env python3
"""Debug the NPMLE: is the fitted mixture actually higher-likelihood than the
truth-implied mixture? And where does the miss mass sit?"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from explore_mosaic import block_rates, load_strata, MOSAIC, BLOCK
from model_mosaic import (V_GRID, A_GRID, kernel_tensor, counts_from_strata,
                          npmle_em, load_ani_row)
from tagcount import read_fasta_single, digest_all


def loglik(counts, K, samp):
    ll = 0.0
    for e, c in counts.items():
        P = K[e] @ samp
        P = np.maximum(P, 1e-300)
        ll += float(c @ np.log(P))
    return ll


def check(dirname, strata_dir, qname, rates, block):
    n = len(read_fasta_single(os.path.join(MOSAIC, "ref.fasta")))
    strata = load_strata(os.path.join(strata_dir, qname + ".strata.tsv"))
    counts = counts_from_strata(strata)
    K = kernel_tensor(strata)

    # truth-implied genome-wide weights over the v grid
    v_true = -np.log(1 - rates)
    w_true, _ = np.histogram(v_true, bins=len(V_GRID) * 0 + V_GRID)
    # histogram needs bin edges; build manually
    edges = np.concatenate([[0], 0.5 * (V_GRID[1:] + V_GRID[:-1]), [np.inf]])
    w_true, _ = np.histogram(v_true, bins=edges)
    w_true = w_true / w_true.sum()

    w_fit = npmle_em(counts, K, iters=20000)
    samp_fit = w_fit  # unconstrained: sample == genome weights in this fit

    # truth restricted to chain sample: approximate by weighting with observed
    # coverage — skip; compare genome truth likelihood directly (upper bound
    # reference only)
    ll_fit = loglik(counts, K, samp_fit)
    ll_true_gw = loglik(counts, K, w_true)
    ani_fit = float(w_fit @ A_GRID) * 100
    ani_true = float(w_true @ A_GRID) * 100

    print(f"\n== {qname}")
    print(f"  truth genome ANI (grid) {ani_true:.2f}")
    print(f"  NPMLE ANI {ani_fit:.2f}   ll_fit {ll_fit:.1f}  ll_true(genomewide) {ll_true_gw:.1f}")
    # where is the mass?
    for lo, hi in [(0, 0.02), (0.02, 0.06), (0.06, 0.12), (0.12, 0.25), (0.25, 1.0)]:
        m = (V_GRID >= lo) & (V_GRID < hi)
        print(f"  v [{lo:.2f},{hi:.2f}): fitted {w_fit[m].sum():.3f}  true {w_true[m].sum():.3f}")
    # observed vs fitted category frequencies
    for e, c in counts.items():
        P = K[e] @ samp_fit
        print(f"  {e:<5} obs {np.round(c / c.sum(), 3)}  fit {np.round(P, 3)}")


if __name__ == "__main__":
    ref = read_fasta_single(os.path.join(MOSAIC, "ref.fasta"))
    n = len(ref)
    nb = (n + BLOCK - 1) // BLOCK
    # gamma_a1.0_ani90 is CASES index 4; bimodal_70core_ani90 index 8
    for i, qname in [(4, "q_ani0.9000__gamma_a1.0_ani90"),
                     (8, "q_ani0.8999__bimodal_70core_ani90")]:
        label = ["gamma_a0.5_ani95", "gamma_a1.0_ani95", "gamma_a2.0_ani95",
                 "gamma_a0.5_ani90", "gamma_a1.0_ani90", "gamma_a1.0_ani98",
                 "bimodal_70core_ani95", "bimodal_50core_ani95",
                 "bimodal_70core_ani90"][i]
        regime = "gamma" if label.startswith("gamma") else "bimodal"
        mean_ani = float(label.split("ani")[1]) / 100
        param = float(label.split("_")[1][1:]) if regime == "gamma" else float(label.split("_")[1][:2]) / 100
        rates = block_rates(regime, mean_ani, param, nb, 90000 + i)
        check(MOSAIC, os.path.join(HERE := os.path.dirname(os.path.abspath(__file__)), "strata_mosaic"),
              qname, rates, BLOCK)
