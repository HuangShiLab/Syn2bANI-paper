#!/usr/bin/env python3
"""Final candidate x dataset summary for MODELS.md. Recomputes the sim fits
(19 mosaic + 12 uniform), adds the BIC comparison gamma vs two-component,
reads the GTDB per-pair results from gtdb_spatial.tsv, and writes
summary_table.tsv.
"""
import csv
import os
import sys
from math import exp, log

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from explore_mosaic import load_strata
from model_mosaic import kernel_tensor, counts_from_strata, npmle_em, \
    load_ani_row, collect_cases
import model_mosaic as mm
from model_discrete import discrete_fit
from het_fit import estimate_het, Stratum

HERE = os.path.dirname(os.path.abspath(__file__))
GF = os.path.join(HERE, "..", "gating_flag")
LRT_THR = 5.0

mm.V_GRID = mm.V_GRID[mm.V_GRID <= 0.450001]
mm.A_GRID = np.exp(-mm.V_GRID)


def fit_all(strata):
    counts = counts_from_strata(strata)
    K = kernel_tensor(strata)
    w = npmle_em(counts, K, iters=5000)
    ani_np = 100 * float(w @ mm.A_GRID)
    ll_np = 0.0
    for e, c in counts.items():
        P = K[e] @ w
        ll_np += float(c @ np.log(np.maximum(P, 1e-300)))
    nll_np = -ll_np
    st = [Stratum(e, s["tag_len"], s["body_len"], list(map(int, s["hist"])),
                  int(s["n_miss"])) for e, s in strata.items()]
    h = estimate_het(st)
    vs, ws, ll_d2 = discrete_fit(counts, K, 2)
    ani_d2 = 100 * float(ws @ np.exp(-vs))
    n_tags = sum(int(c.sum()) for c in counts.values())
    bic_g = 2 * h["nll_het"] + 2 * log(max(n_tags, 2))
    bic_d2 = 2 * (-ll_d2) + 3 * log(max(n_tags, 2))
    lrt = 2 * (h["nll_het"] - nll_np)
    return ani_np, ani_d2, lrt, bic_g, bic_d2, h["nll_het"], nll_np


def gated_of(row):
    gamma = float(row["ani"])
    unif = float(row["ani_uniform"])
    gap = abs(float(row["ani_from_loss"]) - float(row["ani_from_hist"]))
    return gamma, unif if (np.isfinite(gap) and gap > 5.0) else gamma


def summarize(name, truth, ests):
    out = {}
    for k, v in ests.items():
        e = np.array(v) - np.array(truth)
        out[k] = (float(np.abs(e).mean()), float(e.mean()))
    return out


def main():
    results = {}
    bic_rows = []

    # ---- mosaic 19 ----
    truth, ests = [], {"gamma": [], "gated": [], "npmle": [], "D2": [],
                       "spatial": [], "avg": []}
    for c in collect_cases():
        ani_np, ani_d2, lrt, bic_g, bic_d2, nll_g, nll_np = fit_all(c["strata"])
        gamma, gated = gated_of(c["row"])
        spatial = ani_np if lrt > LRT_THR else gated
        truth.append(c["truth"] * 100)
        for k, v in [("gamma", gamma), ("gated", gated), ("npmle", ani_np),
                     ("D2", ani_d2), ("spatial", spatial),
                     ("avg", 0.5 * (gated + ani_np))]:
            ests[k].append(v)
        bic_rows.append((c["name"], bic_g, bic_d2,
                         "D2" if bic_d2 < bic_g else "gamma"))
    results["mosaic19"] = summarize("mosaic", truth, ests)
    n_d2 = sum(1 for r in bic_rows if r[3] == "D2")
    print(f"BIC prefers D2 on {n_d2}/19 mosaic cases")

    # ---- uniform 12 ----
    sd = os.path.join(HERE, "simindel_strata")
    truth_u, ests_u = [], {"gamma": [], "gated": [], "npmle": [], "D2": [],
                           "spatial": [], "avg": []}
    with open("/Users/macstudio/Downloads/Syn2bANI/prototype/simindel/manifest.tsv") as fh:
        next(fh)
        man = {l.split("\t")[0]: float(l.split("\t")[1]) * 100
               for l in fh if l.strip()}
    for name in sorted(man):
        strata = load_strata(os.path.join(sd, name + ".strata.tsv"))
        row = load_ani_row(os.path.join(sd, name + ".ani.tsv"))
        ani_np, ani_d2, lrt, *_ = fit_all(strata)
        gamma, gated = gated_of(row)
        spatial = ani_np if lrt > LRT_THR else gated
        truth_u.append(man[name])
        for k, v in [("gamma", gamma), ("gated", gated), ("npmle", ani_np),
                     ("D2", ani_d2), ("spatial", spatial),
                     ("avg", 0.5 * (gated + ani_np))]:
            ests_u[k].append(v)
    results["uniform12"] = summarize("uniform", truth_u, ests_u)

    # ---- GTDB 2053 (from gtdb_spatial.tsv) ----
    truth_g, ests_g = [], {"gated": [], "npmle": [], "spatial": [], "avg": []}
    bands_g, ests_band = [], {"gated": [], "npmle": [], "spatial": [], "avg": []}
    with open(os.path.join(HERE, "gtdb_spatial.tsv")) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            band, anim, gated, np_, lrt, spatial = p[2], float(p[3]), float(p[4]), \
                float(p[5]), float(p[6]), float(p[7])
            truth_g.append(anim)
            bands_g.append(band)
            for k, v in [("gated", gated), ("npmle", np_), ("spatial", spatial),
                         ("avg", 0.5 * (gated + np_))]:
                ests_g[k].append(v)
    results["gtdb2053"] = summarize("gtdb", truth_g, ests_g)

    # by band
    print("\nGTDB by band:")
    bands = ["0.8-0.85", "0.85-0.9", "0.9-0.95", "0.95-0.99"]
    band_tbl = {}
    for b in bands:
        m = [i for i, x in enumerate(bands_g) if x == b]
        sub = {k: [v[i] for i in m] for k, v in ests_g.items()}
        tr = [truth_g[i] for i in m]
        band_tbl[b] = summarize(b, tr, sub)

    # ---- write summary ----
    with open(os.path.join(HERE, "summary_table.tsv"), "w") as fh:
        fh.write("dataset\testimator\tMAE\tbias\n")
        for ds, tbl in results.items():
            for k, (mae, bias) in tbl.items():
                fh.write(f"{ds}\t{k}\t{mae:.3f}\t{bias:+.3f}\n")
        for b, tbl in band_tbl.items():
            for k, (mae, bias) in tbl.items():
                fh.write(f"gtdb:{b}\t{k}\t{mae:.3f}\t{bias:+.3f}\n")
    with open(os.path.join(HERE, "bic_mosaic.tsv"), "w") as fh:
        fh.write("case\tbic_gamma\tbic_D2\tpreferred\n")
        for r in bic_rows:
            fh.write(f"{r[0]}\t{r[1]:.1f}\t{r[2]:.1f}\t{r[3]}\n")

    hdr = f"{'dataset':<10} {'est':<8} {'MAE':>7} {'bias':>7}"
    print(hdr)
    for ds, tbl in results.items():
        for k, (mae, bias) in tbl.items():
            print(f"{ds:<10} {k:<8} {mae:7.3f} {bias:+7.3f}")
    print("\nGTDB bands:")
    for b, tbl in band_tbl.items():
        row = " ".join(f"{k}:{v[0]:.3f}" for k, v in tbl.items())
        print(f"  {b:<10} {row}")


if __name__ == "__main__":
    main()
