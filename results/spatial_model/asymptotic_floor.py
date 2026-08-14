#!/usr/bin/env python3
"""Asymptotic identifiability floor: build EXPECTED in-chain category counts
from the true per-tag divergences (no sampling noise), then fit each candidate
model. Whatever error remains at infinite data is structural — no count-based
estimator can beat it.

Also separates the two bias components per case:
  family  = error of the best count-based fit on the chained sample
  coverage = error from the chained sample missing the divergent tail
"""
import os
import sys
from math import comb, exp, log

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tagcount import read_fasta_single, digest_all
from explore_mosaic import CASES, block_rates, load_paf, load_strata, MOSAIC, BLOCK
from model_mosaic import V_GRID, A_GRID, kernel_tensor, counts_from_strata, \
    npmle_em, load_ani_row, TOL
from model_spatial import fit_2comp, fit_gamma_cov, pi_cov
from het_fit import estimate_het, Stratum

HERE = os.path.dirname(os.path.abspath(__file__))
EXTRA = os.path.join(MOSAIC, "extra")

ENZ_GEOM = {"BcgI": (32, 26), "AlfI": (32, 26), "AloI": (27, 20), "FalI": (27, 21)}


def expected_counts(v_tags, tag_len, body_len, scale=1e6):
    """Expected category counts for tags with per-site divergences v_tags."""
    a = np.exp(-np.asarray(v_tags))
    P = np.zeros(TOL + 2)
    for m in range(TOL + 1):
        P[m] = np.mean(comb(body_len, m) * (1 - a) ** m * a ** (tag_len - m))
    P[TOL + 1] = 1.0 - P[:TOL + 1].sum()
    return P * scale


def fit_gamma_on_counts(counts, geoms):
    """Run the exact Python port of the Rust het fit on expected counts."""
    strata = []
    for e, c in counts.items():
        tl, bl = geoms[e]
        hist = [int(round(x)) for x in c[:TOL + 1]]
        nm = int(round(c[TOL + 1]))
        strata.append(Stratum(e, tl, bl, hist, nm))
    return estimate_het(strata)


def main():
    ref = read_fasta_single(os.path.join(MOSAIC, "ref.fasta"))
    n = len(ref)

    all_cases = []
    for i, (label, regime, mean_ani, param) in enumerate(CASES):
        all_cases.append((MOSAIC, label, regime, param, mean_ani, BLOCK, 90000 + i))
    mp = os.path.join(EXTRA, "manifest.tsv")
    with open(mp) as fh:
        next(fh)
        for line in fh:
            nm, tr, reg, prm, blk = line.rstrip("\n").split("\t")
            seed = int(nm.split("_")[0][1:])
            mean_ani = float(nm.split("_ani")[1].split("_")[0])
            all_cases.append((EXTRA, nm, reg, float(prm), mean_ani, int(blk), seed))

    print(f"{'case':<34} {'truth':>6} {'g_asy':>6} {'2c_asy':>6} {'np_asy':>6} "
          f"{'chainT':>6} {'AF':>5}")
    floors = {}
    for d, name, regime, param, mean_ani, block, seed in all_cases:
        nb = (n + block - 1) // block
        rates = block_rates(regime if regime in ("gamma", "bimodal") else
                            ("gamma" if "gamma" in name else "bimodal"),
                            mean_ani, param, nb, seed)
        v_blocks = -np.log(1 - rates)
        if d == MOSAIC:
            with open(os.path.join(MOSAIC, "manifest.tsv")) as fh:
                next(fh)
                truth = next(float(l.split("\t")[1]) for l in fh
                             if l.split("\t")[0].endswith("__" + name))
                qname = next(l.split("\t")[0] for l in open(os.path.join(MOSAIC, "manifest.tsv")).readlines()[1:]
                             if l.split("\t")[0].endswith("__" + name))
        else:
            qname = name
            truth = None
            with open(mp) as fh:
                next(fh)
                for l in fh:
                    if l.startswith(name + "\t"):
                        truth = float(l.split("\t")[1])

        qseq = read_fasta_single(os.path.join(d, qname + ".fasta"))
        tags = digest_all(qseq)
        spans = load_paf(os.path.join(d, qname + ".paf"))
        cov = np.zeros(n, dtype=bool)
        for lo, hi in spans:
            cov[lo:hi] = True

        # per-tag true divergence; in-chain subset by extended-span coverage
        counts = {}
        v_all = []
        v_in = []
        for e, pos in tags.items():
            pos = np.array(pos)
            vt = v_blocks[np.minimum(pos // block, nb - 1)]
            v_all.append(vt)
            v_in.append(vt[cov[pos]])
        v_all = np.concatenate(v_all)
        v_in = np.concatenate(v_in)
        truth_pts = truth * 100
        chain_truth = 100 * float(np.mean(np.exp(-v_in)))
        af = float(load_ani_row(os.path.join(d, qname + ".ani.tsv"))["af_query"])

        for e, pos in tags.items():
            pos = np.array(pos)
            vt = v_blocks[np.minimum(pos // block, nb - 1)][cov[pos]]
            counts[e] = expected_counts(vt, *ENZ_GEOM[e])

        # 1) gamma fit on asymptotic in-chain counts
        hg = fit_gamma_on_counts(counts, ENZ_GEOM)
        ani_g = hg["ani"] * 100

        # 2) two-component ML on the same
        K = kernel_tensor({e: dict(tag_len=g[0], body_len=g[1])
                           for e, g in ENZ_GEOM.items()})
        v1, v2, w, _ = fit_2comp(counts, K)
        ani_2c = 100 * (w * exp(-v1) + (1 - w) * exp(-v2))

        # 3) NPMLE on the same
        wnp = npmle_em(counts, K, iters=5000)
        ani_np = 100 * float(wnp @ A_GRID)

        print(f"{name:<34} {truth_pts:6.2f} {ani_g:6.2f} {ani_2c:6.2f} {ani_np:6.2f} "
              f"{chain_truth:6.2f} {af:5.3f}")
        for k, val in [("gamma_asy", ani_g), ("2comp_asy", ani_2c),
                       ("npmle_asy", ani_np), ("chain_truth", chain_truth)]:
            floors.setdefault(k, []).append(val - truth_pts)

    print()
    for k, errs in floors.items():
        errs = np.array(errs)
        print(f"{k:<12} MAE {np.abs(errs).mean():6.3f}  bias {errs.mean():+6.3f}")


if __name__ == "__main__":
    main()
