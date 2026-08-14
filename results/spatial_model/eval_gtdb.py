#!/usr/bin/env python3
"""GTDB-ANIm evaluation of the spatial candidates (2,074 pairs).

Per pair: capped NPMLE fit from strata, NLL, LRT against the gamma fit
(nll_het from the gating-flag refit cache), gated selection at LRT > 5.
Writes per-pair estimates to gtdb_spatial.tsv and prints MAE/bias by band.
"""
import os
import sys
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from explore_mosaic import load_strata
from model_mosaic import kernel_tensor, counts_from_strata, npmle_em
import model_mosaic as mm

HERE = os.path.dirname(os.path.abspath(__file__))
GF = os.path.join(HERE, "..", "gating_flag")
STRATA = os.path.join(GF, "strata_2074")

mm.V_GRID = mm.V_GRID[mm.V_GRID <= 0.450001]
mm.A_GRID = np.exp(-mm.V_GRID)

LRT_THR = 5.0


def fit_pair(fname):
    key = fname[:-4]  # strip .tsv
    q_acc, r_acc = key.split("__")
    try:
        strata = load_strata(os.path.join(STRATA, fname))
        counts = counts_from_strata(strata)
        K = kernel_tensor(strata)
        w = npmle_em(counts, K, iters=5000)
        ani_np = 100 * float(w @ mm.A_GRID)
        ll = 0.0
        for e, c in counts.items():
            P = K[e] @ w
            ll += float(c @ np.log(np.maximum(P, 1e-300)))
        return q_acc, r_acc, ani_np, -ll
    except Exception as e:
        return q_acc, r_acc, float("nan"), float("nan")


def main():
    files = sorted(os.listdir(STRATA))
    with Pool(8) as pool:
        results = pool.map(fit_pair, files, chunksize=16)
    fits = {(q, r): (ani, nll) for q, r, ani, nll in results}

    # gamma NLL cache
    nll_het = {}
    with open(os.path.join(GF, "refit_cache.tsv")) as fh:
        header = next(fh).rstrip("\n").split("\t")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            d = dict(zip(header, p))
            try:
                nll_het[(d["q_acc"], d["r_acc"])] = float(d["nll_het"])
            except (ValueError, KeyError):
                pass

    # joined table
    rows = []
    with open(os.path.join(GF, "gtdb_anim_joined.tsv")) as fh:
        header = next(fh).rstrip("\n").split("\t")
        for line in fh:
            d = dict(zip(header, line.rstrip("\n").split("\t")))
            rows.append(d)

    out = open(os.path.join(HERE, "gtdb_spatial.tsv"), "w")
    out.write("query_asm\tref_asm\tband\tanim_ani\tani_gated\tani_npmle\tlrt_np\tani_spatial\n")
    data = []
    for d in rows:
        key = (d["query_asm"], d["ref_asm"])
        if key not in fits:
            continue
        ani_np, nll_np = fits[key]
        try:
            anim = float(d["anim_ani"])
            gated = float(d["ani"])  # note: 'ani' col = gamma; gated needs gate rule
        except ValueError:
            continue
        # reconstruct the gated estimate exactly as RULES.md: gamma unless
        # |loss-hist| > 5 -> uniform
        gamma = float(d["ani"])
        unif = float(d["ani_uniform"])
        try:
            gap = abs(float(d["ani_from_loss"]) - float(d["ani_from_hist"]))
        except ValueError:
            gap = float("nan")
        gated = unif if (np.isfinite(gap) and gap > 5.0) else gamma
        nll_g = nll_het.get(key, float("nan"))
        lrt = 2 * (nll_g - nll_np) if np.isfinite(nll_g) else float("nan")
        ani_spatial = ani_np if (np.isfinite(lrt) and lrt > LRT_THR) else gated
        out.write(f"{key[0]}\t{key[1]}\t{d['band']}\t{anim}\t{gated}\t{ani_np}"
                  f"\t{lrt}\t{ani_spatial}\n")
        data.append(dict(band=d["band"], anim=anim, gated=gated, np=ani_np,
                         lrt=lrt, spatial=ani_spatial, flag=d["flag"]))
    out.close()

    bands = ["0.8-0.85", "0.85-0.9", "0.9-0.95", "0.95-0.99"]
    print(f"n = {len(data)}")
    print(f"{'estimator':<10} {'all':>7} " + " ".join(f"{b:>10}" for b in bands)
          + f" {'bias':>7}")
    for name, get in [("gated", lambda d: d["gated"]),
                      ("npmle", lambda d: d["np"]),
                      ("spatial", lambda d: d["spatial"])]:
        errs = np.array([get(d) - d["anim"] for d in data])
        line = f"{name:<10} {np.abs(errs).mean():7.3f} "
        for b in bands:
            e = np.array([get(d) - d["anim"] for d in data if d["band"] == b])
            line += f" {np.abs(e).mean():10.3f}"
        print(line + f" {errs.mean():+7.3f}")
    # gate fire rate
    fired = sum(1 for d in data if np.isfinite(d["lrt"]) and d["lrt"] > LRT_THR)
    print(f"LRT>{LRT_THR} fired on {fired}/{len(data)}")


if __name__ == "__main__":
    main()
