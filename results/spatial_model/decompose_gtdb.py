#!/usr/bin/env python3
"""GTDB bias decomposition diagnostics.

1. Correlation of gated-estimator error with internal features (which
   quantities carry the bias information).
2. The AF-mixture family on real data: ANI = AF*E_c + (1-AF)*a_u with a fixed
   saturating a_u — documents why the mosaic-sim coverage correction cannot
   transfer (unchained mass on GTDB is accessory-dominated, and ANIm excludes
   accessory too).
3. npmle residual error vs features: is what remains mechanistic at all?
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GF = os.path.join(HERE, "..", "gating_flag")


def main():
    rows = []
    with open(os.path.join(GF, "gtdb_anim_joined.tsv")) as fh:
        header = next(fh).rstrip("\n").split("\t")
        for line in fh:
            d = dict(zip(header, line.rstrip("\n").split("\t")))
            rows.append(d)

    # spatial estimates from eval_gtdb output
    sp = {}
    with open(os.path.join(HERE, "gtdb_spatial.tsv")) as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            sp[(p[0], p[1])] = dict(npmle=float(p[5]), lrt=float(p[6]),
                                    spatial=float(p[7]))

    data = []
    for d in rows:
        try:
            anim = float(d["anim_ani"])
            gamma = float(d["ani"])
            unif = float(d["ani_uniform"])
            af = float(d["af_query"])
            ret = float(d["retention"])
        except (ValueError, KeyError):
            continue
        gap = abs(float(d["ani_from_loss"]) - float(d["ani_from_hist"]))
        gated = unif if (np.isfinite(gap) and gap > 5.0) else gamma
        key = (d["query_asm"], d["ref_asm"])
        if key not in sp:
            continue
        data.append(dict(anim=anim, gated=gated, af=af, ret=ret, gap=gap,
                         band=d["band"],
                         het_shape=d["het_shape"], n_anchors=float(d["n_anchors"]),
                         n_tags=float(d["n_tags"]), n_chains=float(d["n_chains"]),
                         synteny=float(d["synteny_score"]),
                         bp=float(d["breakpoint_count"]),
                         npmle=sp[key]["npmle"], spatial=sp[key]["spatial"]))

    err = np.array([d["gated"] - d["anim"] for d in data])
    feats = {
        "af_query": np.array([d["af"] for d in data]),
        "retention": np.array([d["ret"] for d in data]),
        "gap_lh": np.array([d["gap"] for d in data]),
        "log_n_anchors_per_tag": np.log10(np.array([d["n_anchors"] for d in data])
                                          / np.array([d["n_tags"] for d in data])),
        "synteny_score": np.array([d["synteny"] for d in data]),
        "bp_per_chain": np.array([d["bp"] for d in data]) / np.maximum([d["n_chains"] for d in data], 1),
    }
    print("correlation of gated error with features (all pairs, n=%d):" % len(data))
    for k, v in feats.items():
        m = np.isfinite(v)
        print(f"  {k:<24} r = {np.corrcoef(v[m], err[m])[0,1]:+.3f}")

    # AF-mixture with fixed saturating a_u (mosaic-calibrated 0.70)
    for a_u in (0.70, 0.80, 0.85):
        est = np.array([d["af"] * d["gated"] + (1 - d["af"]) * a_u * 100
                        for d in data])
        e = est - np.array([d["anim"] for d in data])
        print(f"AF-mix a_u={a_u:.2f}: MAE {np.abs(e).mean():.3f} bias {e.mean():+.3f}")

    # same but only among high-AF pairs (where accessory is small)
    for afmin in (0.5, 0.7):
        sel = [d for d in data if d["af"] >= afmin]
        if not sel:
            continue
        est = np.array([d["af"] * d["gated"] + (1 - d["af"]) * 70 for d in sel])
        e_mix = est - np.array([d["anim"] for d in sel])
        e_gated = np.array([d["gated"] - d["anim"] for d in sel])
        print(f"AF>={afmin} (n={len(sel)}): gated MAE {np.abs(e_gated).mean():.3f} "
              f"vs AF-mix(0.70) MAE {np.abs(e_mix).mean():.3f}")

    # npmle residual vs retention: does the family fix interact with the
    # statistical signal?
    e_np = np.array([d["npmle"] - d["anim"] for d in data])
    ret = feats["retention"]
    for lo, hi in [(0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.01)]:
        m = (ret >= lo) & (ret < hi)
        print(f"retention [{lo},{hi}): n={m.sum():4d}  gated bias {err[m].mean():+6.3f} "
              f" npmle bias {e_np[m].mean():+6.3f}  gated MAE {np.abs(err[m]).mean():6.3f} "
              f" npmle MAE {np.abs(e_np[m]).mean():6.3f}")


if __name__ == "__main__":
    main()
