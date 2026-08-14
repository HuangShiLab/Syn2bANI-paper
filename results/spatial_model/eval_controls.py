#!/usr/bin/env python3
"""Negative controls for the spatial candidates: mid-ANI 15 pairs vs ANIm and
oral/gut same-species pairs vs FastANI. The candidate must do no harm."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from explore_mosaic import load_strata
from model_mosaic import kernel_tensor, counts_from_strata, npmle_em, \
    load_ani_row
import model_mosaic as mm

HERE = os.path.dirname(os.path.abspath(__file__))
GF = os.path.join(HERE, "..", "gating_flag")

mm.V_GRID = mm.V_GRID[mm.V_GRID <= 0.450001]
mm.A_GRID = np.exp(-mm.V_GRID)
LRT_THR = 5.0


def gamma_nll_from_strata(strata):
    from het_fit import estimate_het, Stratum
    st = [Stratum(e, s["tag_len"], s["body_len"], list(map(int, s["hist"])),
                  int(s["n_miss"])) for e, s in strata.items()]
    return estimate_het(st)


def npmle_fit(strata):
    counts = counts_from_strata(strata)
    K = kernel_tensor(strata)
    w = npmle_em(counts, K, iters=5000)
    ani = 100 * float(w @ mm.A_GRID)
    ll = 0.0
    for e, c in counts.items():
        P = K[e] @ w
        ll += float(c @ np.log(np.maximum(P, 1e-300)))
    return ani, -ll


def gated_of(row):
    gamma = float(row["ani"])
    unif = float(row["ani_uniform"])
    gap = abs(float(row["ani_from_loss"]) - float(row["ani_from_hist"]))
    return unif if (np.isfinite(gap) and gap > 5.0) else gamma


def main():
    # ---- mid-ANI 15 ----
    print("== mid-ANI 15 pairs vs ANIm ==")
    import csv
    meta = list(csv.DictReader(open(os.path.join(GF, "midani_joined.tsv")),
                               delimiter="\t"))
    errs = {"gamma": [], "gated": [], "npmle": [], "spatial": []}
    print(f"{'pair':<40} {'truth':>6} {'gamma':>7} {'gated':>7} {'npmle':>7} {'lrt':>7}")
    for i, m in enumerate(meta):
        p = os.path.join(HERE, "midani_strata", f"midani_{i}.strata.tsv")
        row = load_ani_row(os.path.join(HERE, "midani_strata", f"midani_{i}.ani.tsv"))
        strata = load_strata(p)
        ani_np, nll_np = npmle_fit(strata)
        h = gamma_nll_from_strata(strata)
        lrt = 2 * (h["nll_het"] - nll_np)
        truth = float(m["anim_ani"])
        gamma = float(row["ani"])
        gated = gated_of(row)
        spatial = ani_np if lrt > LRT_THR else gated
        tag = f"{m['query']}__{m['reference']}"
        print(f"{tag:<40} {truth:6.2f} {gamma:7.2f} {gated:7.2f} {ani_np:7.2f} {lrt:7.1f}")
        for k, v in [("gamma", gamma), ("gated", gated), ("npmle", ani_np),
                     ("spatial", spatial)]:
            errs[k].append(v - truth)
    for k, e in errs.items():
        e = np.array(e)
        print(f"  {k:<8} MAE {np.abs(e).mean():6.3f}  bias {e.mean():+6.3f}")

    # ---- oral/gut ----
    print("\n== oral/gut same-species (label==high) vs FastANI ==")
    strata_by_pair = {}
    cur = None
    with open(os.path.join(HERE, "oral_gut_strata.tsv")) as fh:
        header = next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            key = (p[0], p[1])
            strata_by_pair.setdefault(key, []).append(p)
    joined = list(csv.DictReader(open(os.path.join(GF, "oral_gut_joined.tsv")),
                                 delimiter="\t"))
    errs = {"gamma": [], "gated": [], "npmle": [], "spatial": []}
    n_used = 0
    for d in joined:
        if d["label"] != "high":
            continue
        key = (d["query"], d["reference"])
        if key not in strata_by_pair:
            continue
        rows = strata_by_pair[key]
        strata = {p[2]: dict(tag_len=int(p[3]), body_len=int(p[4]),
                             n_miss=int(p[5]),
                             hist=[int(x) for x in p[6].split(",")])
                  for p in rows}
        ani_np, nll_np = npmle_fit(strata)
        h = gamma_nll_from_strata(strata)
        lrt = 2 * (h["nll_het"] - nll_np)
        truth = float(d["fastani_ani"])
        gamma = float(d["ani"])
        gated = gated_of(d)
        spatial = ani_np if lrt > LRT_THR else gated
        n_used += 1
        for k, v in [("gamma", gamma), ("gated", gated), ("npmle", ani_np),
                     ("spatial", spatial)]:
            errs[k].append(v - truth)
    print(f"n = {n_used}")
    for k, e in errs.items():
        e = np.array(e)
        print(f"  {k:<8} MAE {np.abs(e).mean():6.3f}  bias {e.mean():+6.3f}")


if __name__ == "__main__":
    main()
