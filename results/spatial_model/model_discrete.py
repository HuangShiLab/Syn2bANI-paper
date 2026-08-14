#!/usr/bin/env python3
"""Discrete-K divergence-mixture ML (K=2,3,4 point masses with free positions
and weights) — the flexible-but-smooth middle ground between gamma and the
70-point NPMLE. Also AF-mixture variants with a fixed saturating unchained
identity, and NPMLE with capped grid. Gated on the 19 mosaic sims.
"""
import os
import sys
from math import comb, exp, log

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model_mosaic as mm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tagcount import read_fasta_single, digest_all
from explore_mosaic import load_strata, MOSAIC
from model_mosaic import V_GRID, A_GRID, kernel_tensor, counts_from_strata, \
    npmle_em, load_ani_row, collect_cases, TOL

HERE = os.path.dirname(os.path.abspath(__file__))


def discrete_fit(counts, K, n_comp, starts=8, rng_seed=0):
    """ML fit of n_comp point masses over the v grid. Returns (vs, ws, nll)."""
    enzymes = list(counts.keys())
    rng = np.random.default_rng(rng_seed)

    def nll(theta):
        vs = np.sort(theta[:n_comp])
        lw = theta[n_comp:]
        ws = np.exp(lw - np.max(lw))
        ws /= ws.sum()
        ll = 0.0
        for e in enzymes:
            Ke = K[e]  # (tol+2, J)
            jidx = [int(np.argmin(np.abs(mm.V_GRID - v))) for v in vs]
            P = ws @ Ke[:, jidx].T
            ll += float(counts[e] @ np.log(np.maximum(P, 1e-300)))
        return -ll

    best = None
    base_starts = [(0.02, 0.3), (0.01, 0.15, 0.4), (0.01, 0.08, 0.25, 0.5)]
    for s in range(starts):
        v0 = np.sort(rng.uniform(0.005, 0.5, n_comp))
        if s < 3 and n_comp <= 3:
            v0 = np.array(base_starts[n_comp - 2][:n_comp])
        w0 = np.full(n_comp, 1.0 / n_comp)
        theta0 = np.concatenate([v0, np.log(w0)])
        r = minimize(nll, theta0, method="Nelder-Mead",
                     options=dict(maxiter=6000, xatol=1e-7, fatol=1e-7))
        if best is None or r.fun < best.fun:
            best = r
    vs = np.sort(best.x[:n_comp])
    lw = best.x[n_comp:]
    ws = np.exp(lw - np.max(lw))
    ws /= ws.sum()
    order = np.argsort([np.argmin(np.abs(V_GRID - v)) for v in vs])
    return vs, ws, -best.fun


def main():
    cases = collect_cases()
    print(f"{'case':<34} {'truth':>6} {'gamma':>6} {'D2':>6} {'D3':>6} {'D4':>6} "
          f"{'NPcap':>6} {'D3+AF':>7}")
    sums = {}
    for c in cases:
        row, strata = c["row"], c["strata"]
        counts = counts_from_strata(strata)
        K = kernel_tensor(strata)
        af = float(row["af_query"])
        truth = c["truth"] * 100

        est = {}
        for k in (2, 3, 4):
            vs, ws, _ = discrete_fit(counts, K, k)
            est[f"D{k}"] = 100 * float(ws @ np.exp(-vs))
            if k == 3:
                est["D3+AF"] = 100 * (af * float(ws @ np.exp(-vs)) + (1 - af) * 0.70)

        # NPMLE with grid capped at v=0.45
        from model_mosaic import V_GRID as VG_FULL
        import model_mosaic as mm
        saved = mm.V_GRID, mm.A_GRID
        mm.V_GRID = V_GRID[V_GRID <= 0.450001]
        mm.A_GRID = np.exp(-mm.V_GRID)
        Kc = kernel_tensor(strata)
        wnp = npmle_em(counts, Kc, iters=5000)
        est["NPcap"] = 100 * float(wnp @ mm.A_GRID)
        mm.V_GRID, mm.A_GRID = saved

        print(f"{c['name']:<34} {truth:6.2f} {float(row['ani']):6.2f} "
              f"{est['D2']:6.2f} {est['D3']:6.2f} {est['D4']:6.2f} "
              f"{est['NPcap']:6.2f} {est['D3+AF']:7.2f}")
        for key, val in est.items():
            sums.setdefault(key, []).append(val - truth)
        sums.setdefault("gamma", []).append(float(row["ani"]) - truth)

    print()
    for key, errs in sums.items():
        errs = np.array(errs)
        print(f"{key:<8} MAE {np.abs(errs).mean():6.3f}  bias {errs.mean():+6.3f}")


if __name__ == "__main__":
    main()
