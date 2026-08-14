#!/usr/bin/env python3
"""Diagnose the mosaic-sim bias: what do the chains see, what do they miss?

For each of the 9 mosaic cases (exact per-block truth regenerated with the
simulator's seeds), this measures:
  - chain coverage of each block vs the block's true rate (the ascertainment),
  - the identity over the chained fraction only (what an unbiased
    chain-restricted estimator could ever report) vs the gamma estimate,
  - genome-wide vs in-chain tag match rates,
  - the gap between in-chain sample and genome truth, per enzyme.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tagcount import read_fasta_single, digest_all

HERE = os.path.dirname(os.path.abspath(__file__))
MOSAIC = os.path.join(HERE, "mosaic")
BLOCK = 5000

CASES = [
    ("gamma_a0.5_ani95", "gamma", 0.95, 0.5),
    ("gamma_a1.0_ani95", "gamma", 0.95, 1.0),
    ("gamma_a2.0_ani95", "gamma", 0.95, 2.0),
    ("gamma_a0.5_ani90", "gamma", 0.90, 0.5),
    ("gamma_a1.0_ani90", "gamma", 0.90, 1.0),
    ("gamma_a1.0_ani98", "gamma", 0.98, 1.0),
    ("bimodal_70core_ani95", "bimodal", 0.95, 0.70),
    ("bimodal_50core_ani95", "bimodal", 0.95, 0.50),
    ("bimodal_70core_ani90", "bimodal", 0.90, 0.70),
]


def block_rates(regime, mean_ani, param, n_blocks, seed):
    rng = np.random.default_rng(seed)
    mean_rate = 1.0 - mean_ani
    if regime == "gamma":
        mult = rng.gamma(shape=param, scale=1.0 / param, size=n_blocks)
    else:
        is_cons = rng.random(n_blocks) < param
        mult = np.where(is_cons, 0.1, 1.0)
    mult = mult / mult.mean()
    return np.clip(mean_rate * mult, 0.0, 0.75)


def load_paf(path):
    spans = []
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            spans.append((int(p[2]), int(p[3])))
    return spans


def load_ani(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        row = fh.readline().rstrip("\n").split("\t")
    d = dict(zip(header, row))
    for k in ("ani", "ani_uniform", "af_query", "af_reference", "retention",
              "ani_from_loss", "ani_from_hist", "ani_gated", "het_shape"):
        d[k] = float(d[k])
    for k in ("n_anchors", "n_chains", "n_tags"):
        d[k] = int(d[k])
    return d


def load_strata(path):
    strata = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            enz, tl, bl, nm, hist = p[2], int(p[3]), int(p[4]), int(p[5]), p[6]
            h = [int(x) for x in hist.split(",")]
            strata[enz] = dict(tag_len=tl, body_len=bl, n_miss=nm, hist=h)
    return strata


def coverage_mask(spans, n):
    cov = np.zeros(n, dtype=bool)
    for lo, hi in spans:
        cov[lo:hi] = True
    return cov


def p_found(a, tag_len, body_len, tol=2):
    """P(tag matches with <= tol mismatches) at per-site identity a (binomial)."""
    from math import comb
    return sum(comb(body_len, m) * (1 - a) ** m * a ** (tag_len - m)
               for m in range(tol + 1))


def main():
    ref = read_fasta_single(os.path.join(MOSAIC, "ref.fasta"))
    n = len(ref)
    n_blocks = (n + BLOCK - 1) // BLOCK

    print(f"{'case':<24} {'truth':>6} {'gamma':>6} {'unif':>6} {'ANI_chain_true':>14} "
          f"{'AF':>5} {'tagcov':>6} {'loss':>6} {'hist':>6} {'gwloss':>7}")
    rows = []
    for i, (label, regime, mean_ani, param) in enumerate(CASES):
        rates = block_rates(regime, mean_ani, param, n_blocks, 90000 + i)
        rate_at = np.repeat(rates, BLOCK)[:n]
        # manifest name carries exact truth
        manifest = {}
        with open(os.path.join(MOSAIC, "manifest.tsv")) as fh:
            next(fh)
            for line in fh:
                nm, tr, reg, prm = line.rstrip("\n").split("\t")
                manifest[nm] = float(tr)
        qname = next(k for k in manifest if k.endswith("__" + label))
        truth = manifest[qname]

        ani = load_ani(os.path.join(MOSAIC, qname + ".ani.tsv"))
        spans = load_paf(os.path.join(MOSAIC, qname + ".paf"))
        strata = load_strata(os.path.join(HERE, "strata_mosaic", qname + ".strata.tsv"))

        cov = coverage_mask(spans, n)
        # truth restricted to the chained fraction (block-size weighted)
        w = np.ones(n)
        ani_chain_true = 1.0 - (rate_at[cov] * w[cov]).sum() / cov.sum()
        af = cov.mean()

        # per-enzyme tags on the query genome
        qseq = read_fasta_single(os.path.join(MOSAIC, qname + ".fasta"))
        tags = digest_all(qseq)
        n_q = {e: len(p) for e, p in tags.items()}

        # tag-level coverage and genome-wide expected match rate
        tagcov_all, gw_q_num, gw_q_den = [], 0.0, 0.0
        per_enz = {}
        for e, pos in tags.items():
            pos = np.array(pos)
            intag = cov[pos]  # approx: PAF spans are the extended chain spans
            st = strata.get(e)
            a_block = 1.0 - rate_at[pos]
            q_gw = p_found_vec(a_block, st["tag_len"], st["body_len"]) if st else None
            gw_q_num += q_gw.sum()
            gw_q_den += len(pos)
            per_enz[e] = dict(n_q=len(pos), n_in=intag.sum(),
                              q_gw=q_gw.mean(),
                              q_in=(q_gw[intag].mean() if intag.any() else np.nan),
                              q_out=(q_gw[~intag].mean() if (~intag).any() else np.nan))
            tagcov_all.append(intag.mean())
        tagcov = float(np.mean(tagcov_all))

        # in-chain observed found rate vs genome-wide expected
        found = sum(sum(s["hist"]) for s in strata.values())
        total = sum(sum(s["hist"]) + s["n_miss"] for s in strata.values())
        gw_q = gw_q_num / gw_q_den
        # genome-wide loss-only ANI: invert p_found on the tag-mean geometry
        k_eff = np.mean([s["tag_len"] for s in strata.values()])
        b_eff = np.mean([s["body_len"] for s in strata.values()])
        gw_loss = invert_p_found(gw_q, k_eff, b_eff)

        print(f"{label:<24} {truth*100:6.2f} {ani['ani']:6.2f} {ani['ani_uniform']:6.2f} "
              f"{ani_chain_true*100:14.2f} {af:5.3f} {tagcov:6.3f} "
              f"{ani['ani_from_loss']:6.2f} {ani['ani_from_hist']:6.2f} {gw_loss:7.2f}")
        rows.append(dict(label=label, truth=truth, ani=ani, af=af,
                         ani_chain_true=ani_chain_true, tagcov=tagcov,
                         per_enz=per_enz, gw_q=gw_q, n_q=n_q,
                         found_rate_in=found / total, strata=strata))

    # detail on the worst case
    worst = max(rows, key=lambda r: abs(r["ani"]["ani"] / 100 - r["truth"]))
    print("\nworst case:", worst["label"])
    print(f"  truth {worst['truth']*100:.2f}  gamma {worst['ani']['ani']:.2f}  "
          f"chained-fraction truth {worst['ani_chain_true']*100:.2f}")
    print(f"  AF {worst['af']:.3f}  tag coverage {worst['tagcov']:.3f}  "
          f"n_tags in chains {worst['ani']['n_tags']}  n_anchors {worst['ani']['n_anchors']}")
    for e, d in worst["per_enz"].items():
        print(f"  {e:<5} N_q {d['n_q']:>5}  in-chain {d['n_in']:>5}  "
              f"E[q] gw {d['q_gw']:.3f}  in {d['q_in']:.3f}  out {d['q_out']:.3f}")
    print(f"  observed in-chain found rate {worst['found_rate_in']:.3f}  "
          f"vs genome-wide E[q] {worst['gw_q']:.3f}")


def p_found_vec(a, tag_len, body_len, tol=2):
    from math import comb
    out = np.zeros_like(a)
    for m in range(tol + 1):
        out += comb(body_len, m) * (1 - a) ** m * a ** (tag_len - m)
    return out


def invert_p_found(q, tag_len, body_len, tol=2):
    lo, hi = 0.5, 0.9999
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if p_found(mid, int(round(tag_len)), int(round(body_len)), tol) > q:
            lo = mid
        else:
            hi = mid
    return 100 * 0.5 * (lo + hi)


if __name__ == "__main__":
    main()
