#!/usr/bin/env python3
"""Exact Python port of the syn2bANI Rust estimators (src/core/mle.rs),
fitted from per-enzyme strata (query tags inside chained regions).

Used to prototype likelihood-exact gating rules (LRT / BIC / boundary
mixture) without rerunning the Rust binary. Verified to reproduce the
shipped `ani`, `ani_uniform`, `retention`, `het_shape` columns to ~1e-6.
"""
import math
from math import lgamma, log, exp

A_LO, A_HI = 0.50, 0.999_999
ALPHA_LO, ALPHA_HI = 0.1, 200.0
ALPHA_UNIFORM = 1.0e6
LRT_CRIT = 3.841


def ln_binom(n, m):
    if m > n:
        return -math.inf
    return lgamma(n + 1) - lgamma(m + 1) - lgamma(n - m + 1)


class Stratum:
    __slots__ = ("enzyme", "tag_len", "body_len", "hist", "n_miss")

    def __init__(self, enzyme, tag_len, body_len, hist, n_miss):
        self.enzyme = enzyme
        self.tag_len = min(tag_len, 32)
        self.body_len = body_len
        self.hist = hist
        self.n_miss = n_miss

    def total(self):
        return sum(self.hist) + self.n_miss

    def n_found(self):
        return sum(self.hist)


def ln_p_found(a, s, m):
    if m > s.body_len or m > s.tag_len:
        return -math.inf
    return ln_binom(s.body_len, m) + m * log(1 - a) + (s.tag_len - m) * log(a)


def nll(a, strata):
    a = min(max(a, 1e-9), 1.0 - 1e-12)
    total = 0.0
    for s in strata:
        p_found = 0.0
        for m in range(len(s.hist)):
            lp = ln_p_found(a, s, m)
            if math.isfinite(lp):
                p_found += exp(lp)
                if s.hist[m] > 0:
                    total -= s.hist[m] * lp
        if s.n_miss > 0:
            total -= s.n_miss * log(max(1.0 - p_found, 1e-300))
    return total


def minimize(f, lo, hi, iters=200, tol=1e-10):
    INV_PHI = 0.618_033_988_749_894_9
    c = hi - INV_PHI * (hi - lo)
    d = lo + INV_PHI * (hi - lo)
    fc, fd = f(c), f(d)
    for _ in range(iters):
        if abs(hi - lo) < tol:
            break
        if fc < fd:
            hi = d; d = c; fd = fc
            c = hi - INV_PHI * (hi - lo); fc = f(c)
        else:
            lo = c; c = d; fc = fd
            d = lo + INV_PHI * (hi - lo); fd = f(d)
    return 0.5 * (lo + hi)


def expected_retention(a, strata):
    tot = sum(s.total() for s in strata)
    if tot <= 0:
        return math.nan
    acc = 0.0
    for s in strata:
        p = sum(exp(ln_p_found(a, s, m)) for m in range(len(s.hist))
                if math.isfinite(ln_p_found(a, s, m)))
        acc += s.total() * p
    return acc / tot


def estimate_uniform(strata):
    strata = [s for s in strata if s.total() > 0]
    if not strata:
        return math.nan, math.nan, math.nan
    ani = minimize(lambda a: nll(a, strata), A_LO, A_HI)
    h = min(1e-5, max(1.0 - ani, 1e-9) / 4.0)
    f0, fp, fm = nll(ani, strata), nll(ani + h, strata), nll(ani - h, strata)
    curv = (fp - 2 * f0 + fm) / (h * h)
    se = math.sqrt(1.0 / curv) if curv > 0 else math.nan
    return ani, se, f0


def ln_p_found_het(d, alpha, s, m):
    if m > s.body_len:
        return -math.inf
    k, b, mf = float(s.tag_len), float(s.body_len), float(m)
    bd = b * d
    if bd <= 0.0 and m > 0:
        return -math.inf
    ln_bd = 0.0 if m == 0 else log(bd)
    return (mf * ln_bd - lgamma(mf + 1) + alpha * log(alpha) - lgamma(alpha)
            + lgamma(alpha + mf) - (alpha + mf) * log(alpha + d * k))


def nll_het(d, alpha, strata):
    if not (d > 0 and alpha > 0) or not math.isfinite(d) or not math.isfinite(alpha):
        return math.inf
    total = 0.0
    for s in strata:
        p_found = 0.0
        for m in range(len(s.hist)):
            lp = ln_p_found_het(d, alpha, s, m)
            if math.isfinite(lp):
                p_found += exp(lp)
                if s.hist[m] > 0:
                    total -= s.hist[m] * lp
        if s.n_miss > 0:
            total -= s.n_miss * log(max(1.0 - p_found, 1e-300))
    return total if math.isfinite(total) else math.inf


def het_ani(d, alpha):
    return (1.0 + d / alpha) ** (-alpha)


def estimate_het(strata):
    """Mirror of mle::estimate_heterogeneous. Returns dict with ani, d, alpha,
    lrt, supported, nll_het, nll_null, n_tags."""
    strata = [s for s in strata if s.total() > 0]
    n_tags = sum(s.total() for s in strata)
    if not strata or n_tags == 0:
        return dict(ani=math.nan, d=math.nan, alpha=math.nan, lrt=math.nan,
                    supported=False, nll_het=math.nan, nll_null=math.nan, n_tags=0,
                    alpha_raw=math.nan, ani_u=math.nan, se_u=math.nan)
    ani_u, se_u, nll_u = estimate_uniform(strata)

    d_lo, d_hi = 1e-6, 1.0
    best = (math.inf, 0.01, 1.0)

    def scan(lo_ln, hi_ln, steps):
        nonlocal best
        for i in range(steps + 1):
            ln_a = lo_ln + (hi_ln - lo_ln) * i / steps
            alpha = exp(ln_a)
            d = minimize(lambda d: nll_het(d, alpha, strata), d_lo, d_hi)
            v = nll_het(d, alpha, strata)
            if v < best[0]:
                best = (v, d, alpha)

    scan(log(ALPHA_LO), log(ALPHA_HI), 48)
    scan(log(best[2]) - 0.35, log(best[2]) + 0.35, 24)
    nll_best, d, alpha = best

    d_null = minimize(lambda d: nll_het(d, ALPHA_UNIFORM, strata), 1e-6, 1.0)
    nll_null = nll_het(d_null, ALPHA_UNIFORM, strata)
    lrt = max(2.0 * (nll_null - nll_best), 0.0)
    supported = lrt > LRT_CRIT and alpha < ALPHA_HI * 0.99

    return dict(
        ani=het_ani(d, alpha) if supported else ani_u,
        d=d, alpha=alpha if supported else math.inf, alpha_raw=alpha,
        lrt=lrt, supported=supported,
        nll_het=nll_best, nll_null=nll_null, n_tags=n_tags,
        ani_u=ani_u, se_u=se_u,
    )


def load_strata_dir(path):
    """Load the per-pair strata dumps into {(q_acc, r_acc): [Stratum]} keyed by
    the accession pair encoded in the filename."""
    import pathlib
    out = {}
    for f in pathlib.Path(path).glob("*.tsv"):
        q_acc, r_acc = f.stem.split("__")
        strata = []
        with open(f) as fh:
            next(fh)
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 7:
                    continue
                _, _, enz, tl, bl, nm, hist = parts[:7]
                strata.append(Stratum(enz, int(tl), int(bl),
                                      [int(x) for x in hist.split(",")], int(nm)))
        out[(q_acc, r_acc)] = strata
    return out
