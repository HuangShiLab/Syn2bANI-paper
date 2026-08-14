#!/usr/bin/env python3
"""Task A — per-pair gating between the gamma (`ani`) and uniform (`ani_uniform`)
estimators, evaluated against ANIm truth on the GTDB 2,074-pair set.

Only per-pair internal quantities from the current feature matrix are used:
retention, het_shape, std_err, n_anchors, n_chains, n_tags, af_*, and the
gamma-uniform discrepancy. Likelihood-exact rules (LRT/BIC) need per-enzyme
strata and are handled in task_a_strata.py once the strata dump arrives.

Reference points (report §3.5, current matrix): always-gamma 2.881,
always-uniform 3.566 (MAE in ANI points vs ANIm, finite pairs).
"""
import pathlib

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent

BANDS = ["0.8-0.85", "0.85-0.9", "0.9-0.95", "0.95-0.99"]


def mae_by_band(df, est_col, truth="anim_ani"):
    err = (df[est_col] - df[truth]).abs()
    out = {"all": err.mean()}
    for b in BANDS:
        m = df["band"] == b
        out[b] = err[m].mean() if m.any() else np.nan
    out["n"] = len(df)
    return out


def gate(df, rule):
    """Return boolean Series: True = use gamma, False = use uniform."""
    kind, par = rule
    if kind == "always_gamma":
        return pd.Series(True, index=df.index)
    if kind == "always_uniform":
        return pd.Series(False, index=df.index)
    if kind == "retention":
        return df["retention"] >= par
    if kind == "retention_shape":
        t, smin = par
        return (df["retention"] >= t) & (df["het_shape"] > smin)
    if kind == "nanchors":
        return df["n_anchors"] >= par
    if kind == "retention_disc":
        # gamma unless retention low OR the models disagree by more than k SE
        t, k = par
        se = df["std_err"].clip(lower=0.05)
        return (df["retention"] >= t) & ((df["ani"] - df["ani_uniform"]).abs() <= k * se)
    raise ValueError(rule)


def apply_gate(df, rule):
    g = gate(df, rule)
    return df["ani"].where(g, df["ani_uniform"])


def main():
    df = pd.read_csv(HERE / "gtdb_anim_joined.tsv", sep="\t")
    df = df[np.isfinite(df["ani"]) & np.isfinite(df["ani_uniform"])].copy()
    df["het_shape"] = pd.to_numeric(df["het_shape"], errors="coerce")
    print(f"finite pairs: {len(df)}")

    rows = {}
    rows["always_gamma"] = mae_by_band(df, "ani")
    rows["always_uniform"] = mae_by_band(df, "ani_uniform")

    # Oracle: per-pair whichever is closer to truth (upper bound on any gate)
    oracle = np.minimum((df["ani"] - df["anim_ani"]).abs(),
                        (df["ani_uniform"] - df["anim_ani"]).abs())
    ob = {"all": oracle.mean(), "n": len(df)}
    for b in BANDS:
        m = df["band"] == b
        ob[b] = oracle[m].mean()
    rows["oracle"] = ob

    # fraction of pairs where gamma is the oracle choice
    g_better = ((df["ani"] - df["anim_ani"]).abs()
                < (df["ani_uniform"] - df["anim_ani"]).abs())
    print(f"oracle picks gamma on {g_better.mean():.1%} of pairs; by band:")
    for b in BANDS:
        m = df["band"] == b
        print(f"  {b}: gamma better {g_better[m].mean():.1%} (n={m.sum()}), "
              f"median retention {df.loc[m, 'retention'].median():.3f}")

    rules = [("always_gamma", None), ("always_uniform", None)]
    for t in [0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70]:
        rules.append(("retention", t))
    for t in [0.35, 0.45, 0.55]:
        for smin in [0.15, 0.3, 0.5]:
            rules.append(("retention_shape", (t, smin)))
    for n in [200, 400, 600, 800, 1000]:
        rules.append(("nanchors", n))
    for t in [0.40, 0.50]:
        for k in [2, 3, 5]:
            rules.append(("retention_disc", (t, k)))

    for r in rules:
        est = apply_gate(df, r)
        tmp = df.assign(gated=est)
        rows[str(r)] = mae_by_band(tmp, "gated")

    tab = pd.DataFrame(rows).T[["all", *BANDS, "n"]]
    tab = tab.sort_values("all")
    pd.set_option("display.width", 160)
    print("\nMAE vs ANIm (ANI points), gated estimator choice:")
    print(tab.round(3).to_string())
    tab.to_csv(HERE / "task_a_gate_sweep.tsv", sep="\t")

    # --- sanity: winner behaviour on the other two datasets ---
    best_simple = None
    for name in tab.index:
        if "retention" in name or "nanchors" in name:
            best_simple = name
            break
    print(f"\nbest simple rule: {best_simple}")

    mid = pd.read_csv(HERE / "midani_joined.tsv", sep="\t")
    mid = mid.rename(columns={"s2b_ani": "ani"})
    mid["het_shape"] = np.nan
    mid["band"] = "mid"
    for name in ["('always_gamma', None)", "('always_uniform', None)", best_simple]:
        rule = eval(name)
        est = apply_gate(mid, rule)
        err = (est - mid["anim_ani"]).abs()
        print(f"mid-ANI {name}: MAE {err.mean():.3f}  bias {(est - mid['anim_ani']).mean():+.3f}"
              f"  (gamma chosen {gate(mid, rule).mean():.0%})")

    og = pd.read_csv(HERE / "oral_gut_joined.tsv", sep="\t")
    og = og[og["label"] == "high"].copy()
    og["het_shape"] = pd.to_numeric(og["het_shape"], errors="coerce")
    og["band"] = "high"
    for name in ["('always_gamma', None)", "('always_uniform', None)", best_simple]:
        rule = eval(name)
        est = apply_gate(og, rule)
        err = (est - og["fastani_ani"]).abs()
        print(f"oral/gut(high) {name}: MAE {err.mean():.3f}  "
              f"bias {(est - og['fastani_ani']).mean():+.3f}  "
              f"(gamma chosen {gate(og, rule).mean():.0%})")


if __name__ == "__main__":
    main()
