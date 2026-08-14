#!/usr/bin/env python3
"""Task B (post-gating) — flag statistic vs the GATED estimator's error.

The Task A gate (uniform iff |ani_from_loss - ani_from_hist| > 5 points) is
applied first; the flag must then rank the *residual* unreliability and,
crucially, transfer: the same threshold must not invert (flagged MAE >= kept
MAE) on GTDB-ANIm, oral/gut (FastANI) and mid-ANI (ANIm).
"""
import pathlib

import numpy as np
import pandas as pd

from task_b_flag import STATS, add_derived, auc

HERE = pathlib.Path(__file__).resolve().parent
GATE_T = 5.0


def load_all():
    g = pd.read_csv(HERE / "gtdb_anim_joined.tsv", sep="\t")
    g = g[np.isfinite(g["ani"])].copy()
    g = add_derived(g)
    g["gated"] = g["ani"].where(g["gap_lh"] <= GATE_T, g["ani_uniform"])
    g["err"] = (g["gated"] - g["anim_ani"]).abs()

    o = pd.read_csv(HERE / "oral_gut_joined.tsv", sep="\t")
    o = o[o["label"] == "high"].copy()
    o = add_derived(o)
    o["gated"] = o["ani"].where(o["gap_lh"] <= GATE_T, o["ani_uniform"])
    o["err"] = (o["gated"] - o["fastani_ani"]).abs()

    m = pd.read_csv(HERE / "midani_joined.tsv", sep="\t").rename(columns={"s2b_ani": "ani"})
    v = pd.read_csv(HERE / "midani_15_verbose.tsv", sep="\t")
    a2s = pd.read_csv(HERE / "midani_acc2seqid.tsv", sep="\t", header=None,
                      names=["acc", "seq"])
    s2a = dict(zip(a2s["seq"], a2s["acc"]))
    v["query_acc"] = v["query"].map(s2a)
    v["reference_acc"] = v["reference"].map(s2a)
    keep = ["query_acc", "reference_acc", "het_shape", "enzyme_spread",
            "enzyme_chi2", "breakpoint_count", "synteny_score",
            "max_block_anchors", "mean_block_anchors"]
    m = m.merge(v[keep], left_on=["query", "reference"],
                right_on=["query_acc", "reference_acc"], how="left")
    m = add_derived(m)
    m["gated"] = m["ani"].where(m["gap_lh"] <= GATE_T, m["ani_uniform"])
    m["err"] = (m["gated"] - m["anim_ani"]).abs()
    return g, o, m


def transfer_table(g, o, m, score_fn, threshold, direction=">"):
    out = {}
    for name, d in [("gtdb", g), ("oral", o), ("mid", m)]:
        s = score_fn(d)
        fl = (s > threshold) if direction == ">" else (s < threshold)
        fl &= np.isfinite(s)
        kept, flagged = d.loc[~fl, "err"], d.loc[fl, "err"]
        inv = (len(flagged) > 0 and flagged.mean() < kept.mean())
        out[name] = (kept.mean(), len(kept),
                     flagged.mean() if len(flagged) else np.nan, len(flagged), inv)
    return out


def main():
    g, o, m = load_all()
    print(f"gated MAE: gtdb {g['err'].mean():.3f}  oral {o['err'].mean():.3f}  "
          f"mid {m['err'].mean():.3f}")
    y = (g["err"] > 1).astype(int)
    print(f"GTDB frac |err_gated|>1: {y.mean():.1%}")

    rows = []
    for s in STATS:
        v = pd.to_numeric(g[s], errors="coerce")
        ok = np.isfinite(v)
        rows.append(dict(stat=s, auc_gated=auc(v[ok], y[ok])))
    tab = pd.DataFrame(rows).sort_values("auc_gated", ascending=False)
    print("\nAUC for |err_gated|>1 on GTDB:")
    print(tab.round(3).to_string(index=False))
    tab.to_csv(HERE / "task_b_auc_gated.tsv", sep="\t", index=False)

    # transfer of the top few at a couple of thresholds
    for s in tab["stat"].head(6):
        v = pd.to_numeric(g[s], errors="coerce")
        for q in [0.75, 0.85]:
            t = v.quantile(q)
            r = transfer_table(g, o, m, lambda d: pd.to_numeric(d[s], errors="coerce"), t)
            inv_any = any(x[4] for x in r.values())
            cells = " | ".join(
                f"{n}: kept {a:.3f}({na}) flagged {b:.3f}({nb}){' INV' if inv else ''}"
                for n, (a, na, b, nb, inv) in r.items())
            print(f"{s:16s} q{int(q*100)} t={t:8.3g} {'INVERTS' if inv_any else 'ok     '} {cells}")


if __name__ == "__main__":
    main()
