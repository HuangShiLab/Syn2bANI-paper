#!/usr/bin/env python3
"""Calibration v6: extend v5 training set with the new high-ANI GTDB-R207 pairs.

Training set:
  - v5 training rows (2,074 post-rescue + 467 hi95 = 2,520 pairs)
  - high-ANI train split from results/gtdb50k/high_ani_results.tsv
    (genome-level 60/40 split, ANI>=95, disjoint from the v5 genomes)

Test set (for reporting; not used in training):
  - 43,334 held-out GTDB-R207 pairs (results/gtdb50k/s2b_50k.tsv + truth_50k.tsv)
  - high-ANI test split

Protocol: same band-holdout CV + RidgeCV as v5, with the same v5a/v5b feature
set selection rule. Final model written to results/panel_by_band/linear_cal_v6.json.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

import calibration_v2 as c2
import calibration_v3 as c3
import calibration_v4 as c4
import calibration_v5 as c5

ROOT = Path(__file__).resolve().parent.parent
HIGH_RES = ROOT / "results/gtdb50k/high_ani_results.tsv"
S2B_50K = ROOT / "results/gtdb50k/s2b_50k.tsv"
TRUTH_50K = ROOT / "results/gtdb50k/truth_50k.tsv"
PAIRS_50K = ROOT / "results/gtdb50k/pairs_50k.tsv"
OUT_CV = ROOT / "results/panel_by_band/calibration_v6_cv.tsv"
OUT_EXT = ROOT / "results/panel_by_band/calibration_v6_external.tsv"
OUT_JSON = ROOT / "results/panel_by_band/linear_cal_v6.json"

TOP_BAND = c4.TOP_BAND
SEEDS = c4.SEEDS
PRIMARY_SEED = c4.PRIMARY_SEED


def load_high_ani_train():
    """Load high-ANI train split and return a matrix compatible with v5."""
    df = pd.read_csv(HIGH_RES, sep="\t")
    df = df[df["split"] == "train"].copy()
    df = df.replace([np.inf, -np.inf], np.nan)
    # recompute gate rule from raw features to guard against column drift
    df = c3.apply_gate_rule(df)
    # band label compatible with eval_pairs
    df["band"] = df["anim_ani"].apply(
        lambda x: TOP_BAND if x >= 95.0 else "other")
    # keep accession identifiers alongside seqids for traceability
    df = df.rename(columns={"q_acc": "query_asm", "r_acc": "ref_asm"})
    df["query_asm"] = df["query_asm"].astype(str)
    df["ref_asm"] = df["ref_asm"].astype(str)
    n_drop = int(df["ani_gated"].isna().sum())
    df = df[df["ani_gated"].notna()].reset_index(drop=True)
    print(f"high-ANI train rows: {len(df)} (dropped {n_drop} non-finite)")
    return df


def load_heldout_50k():
    """Return 43,334 held-out pairs with v5-calibrated and raw gated scores."""
    s2b = pd.read_csv(S2B_50K, sep="\t")
    truth = pd.read_csv(TRUTH_50K, sep="\t")
    pairs = pd.read_csv(PAIRS_50K, sep="\t")
    pairs["pairid"] = pairs["q_acc"] + "__" + pairs["r_acc"]
    df = truth.merge(pairs[["pairid", "band"]], on="pairid", how="inner")
    df = df.merge(s2b, on="pairid", how="inner")
    df = df.replace([np.inf, -np.inf], np.nan)
    df["gate_fallback"] = (df["gate"] == "uniform_fallback").astype(float)
    return df


def evaluate_on(df, fitted, features, label):
    """Compute per-band and overall MAE/bias/r for a test set."""
    df = df.copy()
    df["cal"] = c2.predict_dataframe(fitted, df, features)
    rows = []
    for method, col in [("raw_gated", "ani_gated"), ("calibrated_v6", "cal")]:
        if col not in df.columns or df[col].isna().all():
            continue
        sub = df[df[col].notna() & df["anim_ani"].notna()].copy()
        for band in sorted(sub["band"].unique()):
            b = sub[sub["band"] == band]
            err = b[col] - b["anim_ani"]
            rows.append({
                "dataset": label, "band": band, "method": method,
                "n": len(b), "MAE": float(err.abs().mean()),
                "bias": float(err.mean()),
                "r": float(b[col].corr(b["anim_ani"])),
            })
        err = sub[col] - sub["anim_ani"]
        rows.append({
            "dataset": label, "band": "all", "method": method,
            "n": len(sub), "MAE": float(err.abs().mean()),
            "bias": float(err.mean()),
            "r": float(sub[col].corr(sub["anim_ani"])),
        })
    return pd.DataFrame(rows)


def main():
    old = c5.load_gated_v5()
    new = c5.load_hi95_v5()
    c5.verify_gate_rule(old, "2074 v9rescue")
    c5.verify_gate_rule(new, "hi95 v9rescue")
    v5_train = pd.concat([old, new], ignore_index=True)
    high_train = load_high_ani_train()

    # ensure compatible columns
    for col in c3.V3B:
        if col not in v5_train.columns:
            v5_train[col] = np.nan
        if col not in high_train.columns:
            high_train[col] = np.nan

    df = pd.concat([v5_train, high_train], ignore_index=True)
    print(f"v6 merged training rows: {len(df)} (v5 {len(v5_train)} + high {len(high_train)})")
    print(f"top band rows: {int((df['band'] == TOP_BAND).sum())}")

    rows = []
    d0 = c4.shuffled(df, PRIMARY_SEED)
    for r in c2.metric_rows("raw_gated", d0["ani_gated"].values,
                            d0["anim_ani"].values, d0["band"].values):
        r["seed"] = PRIMARY_SEED
        rows.append(r)

    # v5 reference on the old training set
    rows += c4.cv_experiment(v5_train, "v5b_ref", c3.V3B, SEEDS)
    # v6 variants on the merged set
    rows += c4.cv_experiment(df, "v6a_gated9", c3.V3A, SEEDS)
    rows += c4.cv_experiment(df, "v6b_gated9+gate", c3.V3B, SEEDS)

    cv = pd.DataFrame(rows)
    cv = c4.add_mean_rows(cv).round(4)
    cv = cv[["experiment", "seed", "band", "n", "MAE", "RMSE", "bias", "r"]]
    cv.to_csv(OUT_CV, sep="\t", index=False)
    pd.set_option("display.width", 160)
    print(cv[cv.band == "all"].to_string(index=False))

    mae_a = c4.overall(cv, "v6a_gated9")
    mae_b = c4.overall(cv, "v6b_gated9+gate")
    use_b = mae_b < 0.995 * mae_a
    feats = c3.V3B if use_b else c3.V3A
    names = c3.RUST_NAMES_V3B if use_b else c3.RUST_NAMES_V3B[:9]
    tag = "v6b" if use_b else "v6a"
    print(f"\nv6a={mae_a:.4f} v6b={mae_b:.4f} (seed {PRIMARY_SEED}) -> final: {tag}")

    cal, fitted, mae_in, alpha = c4.train_final(df, feats, names, PRIMARY_SEED)
    cal["name"] = "gtdb_r207_linear_cal_v6"
    with open(OUT_JSON, "w") as fh:
        json.dump(cal, fh, indent=2)
    print(f"final v6 in-sample MAE {mae_in:.4f}, alpha {alpha}; wrote {OUT_JSON}")

    # external / held-out evaluations
    held = load_heldout_50k()
    high_test = pd.read_csv(HIGH_RES, sep="\t")
    high_test = high_test[high_test["split"] == "test"].copy()
    high_test = high_test.replace([np.inf, -np.inf], np.nan)
    high_test = c3.apply_gate_rule(high_test)
    high_test["band"] = high_test["anim_ani"].apply(
        lambda x: TOP_BAND if x >= 95.0 else "other")

    ext_rows = []
    ext_rows.append(evaluate_on(held, fitted, feats, "gtdb_43k_heldout"))
    ext_rows.append(evaluate_on(high_test, fitted, feats, "high_ani_test"))
    # keep oral/gut and mid-ANI external gates from v5 for continuity
    og_rows, _ = c5.external_oralgut_v5(fitted, feats, "calibrated_v6")
    mid_rows, _ = c5.external_midani_v5(fitted, feats, "calibrated_v6")
    ext = pd.concat(ext_rows + [pd.DataFrame(og_rows), pd.DataFrame(mid_rows)],
                    ignore_index=True).round(4)
    ext.to_csv(OUT_EXT, sep="\t", index=False)
    print("\n=== external validation (v6) ===")
    print(ext.to_string(index=False))


if __name__ == "__main__":
    main()
