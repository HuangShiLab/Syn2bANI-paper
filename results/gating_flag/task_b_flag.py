#!/usr/bin/env python3
"""Task B — which per-pair statistic best predicts large |error| vs truth?

Datasets (truth in ANI points):
  GTDB-ANIm   gtdb_anim_joined.tsv   truth = anim_ani      (primary; AUC here)
  oral/gut    oral_gut_joined.tsv    truth = fastani_ani   (label == high only)
  mid-ANI     midani_joined.tsv      truth = anim_ani      (15 pairs, 87.6-90.2)

For every candidate statistic: AUC for |err| > 1 on GTDB, then the behaviour
at a transfer threshold (kept vs flagged MAE on all three datasets). A flag
"inverts" on a dataset if flagged MAE < kept MAE — the current flag does
exactly that on GTDB (4.11 ok vs 1.97 INCONSISTENT).
"""
import pathlib

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent


def auc(score, label):
    """Mann-Whitney AUC of score for binary label (handles ties)."""
    s = pd.Series(score).rank()
    pos = s[label == 1]
    n_pos, n_neg = (label == 1).sum(), (label == 0).sum()
    if n_pos == 0 or n_neg == 0:
        return np.nan
    return (pos.sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)


def add_derived(df, ani_col="ani"):
    d = df.copy()
    d["het_shape"] = pd.to_numeric(d.get("het_shape"), errors="coerce")
    d["gap_lh"] = (d["ani_from_loss"] - d["ani_from_hist"]).abs()
    d["gap_lh_se"] = d["gap_lh"] / d["std_err"].clip(lower=1e-6)
    d["gap_gu"] = (d[ani_col] - d["ani_uniform"]).abs()
    d["gap_gu_se"] = d["gap_gu"] / d["std_err"].clip(lower=1e-6)
    d["inv_shape"] = 1.0 / d["het_shape"].replace([np.inf], np.nan)
    d["bp_per_anchor"] = d["breakpoint_count"] / d["n_anchors"].clip(lower=1)
    d["anchors_per_tag"] = d["n_anchors"] / d["n_tags"].clip(lower=1)
    d["af_min"] = d[["af_query", "af_reference"]].min(axis=1)
    d["gap_lh_af"] = d["gap_lh"] / d["af_min"].clip(lower=0.01)
    return d


STATS = ["enzyme_chi2", "enzyme_spread", "anchor_adjacency", "breakpoint_count",
         "std_err", "retention", "af_query", "af_reference", "het_shape",
         "inv_shape", "n_anchors", "n_chains", "gap_lh", "gap_lh_se",
         "gap_gu", "gap_gu_se", "bp_per_anchor", "anchors_per_tag", "af_min",
         "gap_lh_af"]


def load_all():
    g = pd.read_csv(HERE / "gtdb_anim_joined.tsv", sep="\t")
    g = g[np.isfinite(g["ani"])].copy()
    g = add_derived(g)
    g["err"] = (g["ani"] - g["anim_ani"]).abs()

    o = pd.read_csv(HERE / "oral_gut_joined.tsv", sep="\t")
    o = o[o["label"] == "high"].copy()
    o = add_derived(o)
    o["err"] = (o["ani"] - o["fastani_ani"]).abs()

    m = pd.read_csv(HERE / "midani_joined.tsv", sep="\t")
    m = m.rename(columns={"s2b_ani": "ani"})
    # merge the verbose 4e features (het_shape, enzyme stats) via seqid map
    v = pd.read_csv(HERE / "midani_15_verbose.tsv", sep="\t")
    a2s = pd.read_csv(HERE / "midani_acc2seqid.tsv", sep="\t", header=None,
                      names=["acc", "seq"])
    s2a = dict(zip(a2s["seq"], a2s["acc"]))
    v["query_acc"] = v["query"].map(s2a)
    v["reference_acc"] = v["reference"].map(s2a)
    keep = ["query_acc", "reference_acc", "het_shape", "enzyme_spread",
            "enzyme_chi2", "breakpoint_count", "anchor_adjacency",
            "max_block_anchors", "mean_block_anchors"]
    m = m.merge(v[keep], left_on=["query", "reference"],
                right_on=["query_acc", "reference_acc"], how="left")
    m = add_derived(m)
    m["err"] = (m["ani"] - m["anim_ani"]).abs()
    return g, o, m


def main():
    g, o, m = load_all()
    print(f"GTDB n={len(g)} (|err|>1: {(g['err'] > 1).mean():.1%}), "
          f"oral/gut n={len(o)}, mid-ANI n={len(m)}")

    y = (g["err"] > 1).astype(int)
    rows = []
    for s in STATS:
        v = pd.to_numeric(g[s], errors="coerce")
        ok = np.isfinite(v)
        rows.append(dict(stat=s, auc=auc(v[ok], y[ok]), frac_finite=ok.mean()))
    tab = pd.DataFrame(rows).sort_values("auc", ascending=False)
    print("\nAUC for |err|>1 on GTDB-ANIm (higher = flags the bad pairs):")
    print(tab.round(3).to_string(index=False))
    tab.to_csv(HERE / "task_b_auc.tsv", sep="\t", index=False)

    # transfer check for the top statistics: flag = stat > t, t = GTDB 85th pct
    print("\nTransfer at GTDB 85th-percentile threshold "
          "(flagged MAE must be >= kept MAE everywhere):")
    for s in tab["stat"].head(10):
        gv = pd.to_numeric(g[s], errors="coerce")
        t = gv.quantile(0.85)
        line = [f"{s:16s} t={t:8.3g}"]
        for name, d in [("gtdb", g), ("oral", o), ("mid", m)]:
            v = pd.to_numeric(d[s], errors="coerce")
            fl = v > t
            fl &= np.isfinite(v)
            kept = d.loc[~fl, "err"]; flagged = d.loc[fl, "err"]
            inv = "INVERTS" if len(flagged) and flagged.mean() < kept.mean() else ""
            line.append(f"{name}: kept {kept.mean():.3f}(n={len(kept)}) "
                        f"flagged {flagged.mean() if len(flagged) else float('nan'):.3f}"
                        f"(n={len(flagged)}) {inv}")
        print("  " + " | ".join(line))


if __name__ == "__main__":
    main()
