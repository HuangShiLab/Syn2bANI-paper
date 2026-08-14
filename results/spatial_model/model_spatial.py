#!/usr/bin/env python3
"""Candidate spatial/ascertainment models, gated on the 19 exact-truth mosaic
sims. All fits use only per-pair strata + feature row (plus exact N_q from
Python tag extraction, which the Rust binary does not yet report).

Candidates:
  C1   two-component ML on in-chain counts (no spatial terms)
  C2   two-component ML + coverage-weighted likelihood + AF constraint
  C3   gamma (d, alpha) + coverage-weighted likelihood + AF constraint
  A1   AF-mixture in identity space: AF*E_c[a] + (1-AF)*a_u, a_u from the
       anchor-residual upper bound
  A2   same in divergence space
  Baselines: gamma (ani), uniform, gated.

pi(v) = 1/(1+exp((v - v50)/s)): the chain-coverage function, an algorithm
property. Frozen constants, sensitivity-tested.
"""
import os
import sys
from math import comb, log, exp, lgamma

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tagcount import read_fasta_single, digest_all
from explore_mosaic import load_strata, MOSAIC
from model_mosaic import V_GRID, A_GRID, kernel_tensor, counts_from_strata, \
    npmle_em, load_ani_row, collect_cases, TOL

HERE = os.path.dirname(os.path.abspath(__file__))

# coverage function constants (algorithm property, calibrated on sims only)
PI_V50 = 0.20
PI_S = 0.06
# multi-match inflation of n_anchors (repeat/duplicate matches), sim-calibrated
KAPPA = 1.16


def pi_cov(v, v50=PI_V50, s=PI_S):
    return 1.0 / (1.0 + np.exp((v - v50) / s))


def q_found_vec(a, tag_len, body_len, tol=TOL):
    out = np.zeros_like(a, dtype=float)
    for m in range(tol + 1):
        out += comb(body_len, m) * (1 - a) ** m * a ** (tag_len - m)
    return out


# ── two-component model ──────────────────────────────────────────────────────

def cat_probs_2comp(v1, v2, w, K, pi=None):
    """Category probabilities per enzyme for a 2-component divergence mixture.
    If pi given, weights are coverage-tilted (in-chain sample)."""
    a1, a2 = exp(-v1), exp(-v2)
    j1 = int(np.argmin(np.abs(V_GRID - v1)))
    j2 = int(np.argmin(np.abs(V_GRID - v2)))
    w1, w2 = w, 1.0 - w
    if pi is not None:
        w1 *= pi[j1]
        w2 *= pi[j2]
    z = w1 + w2
    w1, w2 = w1 / z, w2 / z
    return {e: w1 * K[e][:, j1] + w2 * K[e][:, j2] for e in K}


def nll_2comp(params, counts, K, pi=None, af=None, af_n=0.0):
    lv1, lv2, lw = params
    v1, v2 = np.exp(lv1), np.exp(lv2)
    if v2 < v1:
        v1, v2 = v2, v1
    w = 1.0 / (1.0 + np.exp(-lw))
    P = cat_probs_2comp(v1, v2, w, K, pi)
    ll = 0.0
    for e, c in counts.items():
        ll += float(c @ np.log(np.maximum(P[e], 1e-300)))
    if pi is not None and af is not None and af_n > 0:
        cov = w * pi[int(np.argmin(np.abs(V_GRID - v1)))] + \
              (1 - w) * pi[int(np.argmin(np.abs(V_GRID - v2)))]
        cov = min(max(cov, 1e-6), 1 - 1e-6)
        ll += af_n * (af * log(cov) + (1 - af) * log(1 - cov))
    return -ll


def fit_2comp(counts, K, pi=None, af=None, af_n=0.0):
    best = None
    for v10, v20, w0 in [(0.01, 0.2, 0.7), (0.03, 0.3, 0.5), (0.005, 0.1, 0.8),
                         (0.05, 0.05, 0.5), (0.02, 0.5, 0.6)]:
        r = minimize(nll_2comp, [log(v10), log(v20), log(w0 / (1 - w0))],
                     args=(counts, K, pi, af, af_n), method="Nelder-Mead",
                     options=dict(maxiter=4000, xatol=1e-6, fatol=1e-6))
        if best is None or r.fun < best.fun:
            best = r
    lv1, lv2, lw = best.x
    v1, v2 = sorted([exp(lv1), exp(lv2)])
    w = 1.0 / (1.0 + exp(-lw))
    if exp(lv1) > exp(lv2):
        w = 1 - w
    return v1, v2, w, -best.fun


# ── gamma + coverage + AF ────────────────────────────────────────────────────

def ln_nb(d, alpha, tag_len, body_len, m):
    if m > body_len:
        return -np.inf
    k, b, mf = float(tag_len), float(body_len), float(m)
    bd = b * d
    if bd <= 0.0 and m > 0:
        return -np.inf
    ln_bd = 0.0 if m == 0 else log(bd)
    return (mf * ln_bd - lgamma(mf + 1) + alpha * log(alpha) - lgamma(alpha)
            + lgamma(alpha + mf) - (alpha + mf) * log(alpha + d * k))


def nll_gamma_cov(params, counts, strata, pi, af, af_n):
    ld, la = params
    d, alpha = exp(ld), exp(la)
    # coverage-tilted gamma: sample density ∝ g(r) pi(r d)
    # integrate over the v grid: r d = v, gamma density in v
    v = V_GRID[1:]
    # gamma pdf of v with mean d, shape alpha
    logpdf = (alpha * log(alpha / d) + (alpha - 1) * np.log(v)
              - alpha * v / d - lgamma(alpha))
    gv = np.exp(logpdf)
    gv /= gv.sum()
    pv = pi_cov(v)
    samp = gv * pv
    cov = samp.sum()
    samp /= cov
    ll = 0.0
    for e, c in counts.items():
        s = strata[e]
        # category probs under the coverage-tilted mixture: binomial kernel at
        # per-site identity a = exp(-v), tag sharing its region's v
        a = np.exp(-v)
        Pm = np.zeros(TOL + 2)
        for m in range(TOL + 1):
            Pm[m] = float(samp @ (comb(s["body_len"], m) * (1 - a) ** m
                                  * a ** (s["tag_len"] - m)))
        Pm[TOL + 1] = 1.0 - Pm[:TOL + 1].sum()
        ll += float(c @ np.log(np.maximum(Pm, 1e-300)))
    if af is not None and af_n > 0:
        covc = min(max(cov, 1e-6), 1 - 1e-6)
        ll += af_n * (af * log(covc) + (1 - af) * log(1 - covc))
    return -ll


def fit_gamma_cov(counts, strata, pi=True, af=None, af_n=0.0):
    best = None
    for d0, a0 in [(0.05, 1.0), (0.1, 0.5), (0.02, 2.0), (0.15, 1.0)]:
        r = minimize(nll_gamma_cov, [log(d0), log(a0)],
                     args=(counts, strata, pi, af, af_n), method="Nelder-Mead",
                     options=dict(maxiter=2000))
        if best is None or r.fun < best.fun:
            best = r
    d, alpha = exp(best.x[0]), exp(best.x[1])
    ani = (1 + d / alpha) ** (-alpha)
    return d, alpha, ani, -best.fun


# ── closed-form AF mixtures ──────────────────────────────────────────────────

def invert_p_found(q, tag_len, body_len, tol=TOL):
    if q <= 0:
        return 0.5
    if q >= 1:
        return 1.0
    lo, hi = 0.5, 0.99999
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if q_found_vec(np.array([mid]), tag_len, body_len, tol)[0] > q:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    cases = collect_cases()
    # exact N_q per query genome (Python tag extraction)
    nq_cache = {}
    for c in cases:
        qseq = read_fasta_single(os.path.join(c["dir"], c["name"] + ".fasta"))
        nq_cache[c["name"]] = digest_all(qseq)

    print(f"{'case':<34} {'truth':>6} {'gamma':>6} {'C1_2c':>6} "
          f"{'C2_2cAF':>7} {'C3_gAF':>7} {'A1':>6}")
    sums = {}
    for c in cases:
        row, strata = c["row"], c["strata"]
        counts = counts_from_strata(strata)
        K = kernel_tensor(strata)
        af = float(row["af_query"])
        truth = c["truth"] * 100

        # C1: two-component, no spatial
        v1, v2, w, ll1 = fit_2comp(counts, K)
        ani_c1 = 100 * (w * exp(-v1) + (1 - w) * exp(-v2))

        # C2: two-component + coverage + AF
        pi = pi_cov(V_GRID)
        v1b, v2b, wb, ll2 = fit_2comp(counts, K, pi=pi, af=af, af_n=2000.0)
        ani_c2 = 100 * (wb * exp(-v1b) + (1 - wb) * exp(-v2b))

        # C3: gamma + coverage + AF
        d3, al3, ani_c3f, ll3 = fit_gamma_cov(counts, strata, pi=True,
                                              af=af, af_n=2000.0)
        ani_c3 = 100 * ani_c3f

        # A1: AF * ani_gated + (1-AF) * a_u_bound from anchor residual
        n_q = sum(len(p) for p in nq_cache[c["name"]].values())
        n_anchors = int(row["n_anchors"])
        found_in = sum(sum(s["hist"]) for s in strata.values())
        n_in = sum(sum(s["hist"]) + s["n_miss"] for s in strata.values())
        n_u = max(n_q - n_in, 1)
        anchors_out = max(n_anchors / KAPPA - found_in, 0.0)
        q_u = min(anchors_out / n_u, 1.0)
        # tag-geometry-weighted inversion
        a_u = invert_p_found(max(q_u, 1e-6), 30, 23)
        ani_c = float(row["ani_gated"]) / 100
        ani_a1 = 100 * (af * ani_c + (1 - af) * a_u)

        print(f"{c['name']:<34} {truth:6.2f} {float(row['ani']):6.2f} "
              f"{ani_c1:6.2f} {ani_c2:7.2f} {ani_c3:7.2f} {ani_a1:6.2f}")
        for key, val in [("gamma", float(row["ani"])), ("gated", float(row["ani_gated"])),
                         ("C1", ani_c1), ("C2", ani_c2), ("C3", ani_c3), ("A1", ani_a1)]:
            sums.setdefault(key, []).append(val - truth)

    print()
    for key, errs in sums.items():
        errs = np.array(errs)
        print(f"{key:<8} MAE {np.abs(errs).mean():6.3f}  bias {errs.mean():+6.3f}")


if __name__ == "__main__":
    main()
