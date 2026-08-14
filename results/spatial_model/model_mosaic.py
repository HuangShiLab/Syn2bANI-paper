#!/usr/bin/env python3
"""Candidate mechanistic models for the mosaic bias, fitted per pair from
strata + features. Gated on the 19 exact-truth mosaic sims (9 main + 10 extra).

Candidates:
  NPMLE    unconstrained grid NPMLE on in-chain counts (isolates family
           misspecification from coverage ascertainment)
  T1       two-region moment estimator: in-chain NPMLE for the chained mass,
           unchained rate from the genome-wide anchor residual, AF for the
           mass split. No coverage function pi needed.
  S1       grid NPMLE with coverage weighting pi(v) + AF constraint
  A1/A2    closed-form AF-weighted mixtures (identity / divergence space)
"""
import os
import sys
from math import comb, log, exp

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tagcount import read_fasta_single, digest_all
from explore_mosaic import CASES, block_rates, load_paf, load_strata, MOSAIC, BLOCK

HERE = os.path.dirname(os.path.abspath(__file__))
EXTRA = os.path.join(MOSAIC, "extra")

TOL = 2
ENZ_ORDER = ["AlfI", "AloI", "BcgI", "FalI"]

# divergence grid: per-site identity a = exp(-v)
V_GRID = np.concatenate([[0.0], np.logspace(np.log10(1e-4), np.log10(0.8), 70)])
A_GRID = np.exp(-V_GRID)


def q_kernel(a, tag_len, body_len, tol=TOL):
    """P(m=0..tol) and miss for one stratum geometry at identity a (vector)."""
    a = np.asarray(a)
    out = []
    for m in range(tol + 1):
        out.append(comb(body_len, m) * (1 - a) ** m * a ** (tag_len - m))
    found = sum(out)
    out.append(1.0 - found)
    return np.array(out)  # shape (tol+2, len(a))


# kernel tensor per enzyme: K[e][category, j]
def kernel_tensor(strata):
    K = {}
    for e, s in strata.items():
        K[e] = q_kernel(A_GRID, s["tag_len"], s["body_len"])
    return K


def npmle_em(counts, K, iters=2000, w0=None, pi=None, af=None, af_strength=0.0,
             verbose=False):
    """EM for the grid mixture. counts[e] = category counts (len tol+2).
    If pi is given, the in-chain kernel is weighted: tag sample weights are
    w*pi/sum(w*pi); af constraint adds a binomial pseudo-observation on
    sum(w*pi) = af with af_strength pseudo-counts.
    Returns w (genome-wide weights)."""
    J = len(A_GRID)
    w = np.full(J, 1.0 / J) if w0 is None else w0.copy()
    enzymes = list(counts.keys())
    total_counts = sum(c.sum() for c in counts.values())
    for it in range(iters):
        # E-step responsibilities for in-chain categories
        if pi is not None:
            cov = (w * pi).sum()
            samp = w * pi / max(cov, 1e-300)
        else:
            samp = w
        num = np.zeros(J)
        for e in enzymes:
            P = K[e] @ samp            # category probs, len tol+2
            P = np.maximum(P, 1e-300)
            c = counts[e]
            # responsibility mass per grid point: samp_j * sum_m K_mj c_m / P_m
            num += (samp[None, :] * K[e]).T @ (c / P)
        w_new_samp = num / total_counts  # target for the *sample* weights
        if pi is None:
            w_new = w_new_samp
        else:
            # invert sample -> genome weights approximately: w ∝ samp/pi,
            # then one gradient step on the AF binomial pseudo-likelihood
            w_new = w_new_samp / np.maximum(pi, 1e-12)
            w_new /= w_new.sum()
            if af is not None and af_strength > 0:
                cov = (w_new * pi).sum()
                grad = af_strength * (af - cov) / max(cov * (1 - cov), 1e-6) * pi
                w_new = w_new * np.exp(grad / max(af_strength, 1.0))
                w_new = np.maximum(w_new, 0)
                w_new /= w_new.sum()
        delta = np.abs(w_new - w).max()
        w = w_new
        if delta < 1e-10:
            break
    return w


def counts_from_strata(strata):
    return {e: np.array(s["hist"] + [s["n_miss"]], dtype=float)
            for e, s in strata.items()}


def ani_of(w):
    return float(w @ A_GRID)


def load_ani_row(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        row = dict(zip(header, fh.readline().rstrip("\n").split("\t")))
    return row


def collect_cases():
    """Return list of dicts: name, truth, ani row, strata, N_q (exact)."""
    cases = []
    with open(os.path.join(MOSAIC, "manifest.tsv")) as fh:
        next(fh)
        for line in fh:
            nm, tr, reg, prm = line.rstrip("\n").split("\t")
            cases.append((MOSAIC, nm, float(tr)))
    mp = os.path.join(EXTRA, "manifest.tsv")
    if os.path.exists(mp):
        with open(mp) as fh:
            next(fh)
            for line in fh:
                nm, tr, reg, prm, blk = line.rstrip("\n").split("\t")
                cases.append((EXTRA, nm, float(tr)))
    out = []
    for d, nm, tr in cases:
        strata = load_strata(os.path.join(d if d == EXTRA else HERE + "/strata_mosaic",
                                          nm + ".strata.tsv")
                             if d == EXTRA else os.path.join(HERE, "strata_mosaic", nm + ".strata.tsv"))
        row = load_ani_row(os.path.join(d, nm + ".ani.tsv"))
        out.append(dict(dir=d, name=nm, truth=tr, row=row, strata=strata))
    return out


def main():
    cases = collect_cases()
    print(f"{len(cases)} cases")
    hdr = (f"{'case':<34} {'truth':>6} {'gamma':>6} {'unif':>6} {'gated':>6} "
           f"{'NPMLE':>6} {'T1':>6}")
    print(hdr)
    results = []
    for c in cases:
        row, strata = c["row"], c["strata"]
        counts = counts_from_strata(strata)
        K = kernel_tensor(strata)
        w = npmle_em(counts, K)
        ani_np = ani_of(w) * 100
        results.append(dict(c, ani_npmle=ani_np))
        print(f"{c['name']:<34} {c['truth']*100:6.2f} {float(row['ani']):6.2f} "
              f"{float(row['ani_uniform']):6.2f} {float(row['ani_gated']):6.2f} "
              f"{ani_np:6.2f}")
    # summary
    for key, get in [("gamma", lambda r: float(r["row"]["ani"])),
                     ("uniform", lambda r: float(r["row"]["ani_uniform"])),
                     ("gated", lambda r: float(r["row"]["ani_gated"])),
                     ("npmle", lambda r: r["ani_npmle"])]:
        errs = np.array([get(r) - r["truth"] * 100 for r in results])
        print(f"{key:<8} MAE {np.abs(errs).mean():6.3f}  bias {errs.mean():+6.3f}")


if __name__ == "__main__":
    main()
