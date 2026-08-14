#!/usr/bin/env python3
"""Calibration v3: ridge recalibrated against the gated estimator
(Syn2bANI @ 98177dc; `ani_gated` = gamma unless |ani_from_loss -
ani_from_hist| > 5 points -> uniform fallback; rule verified 2074/2074
against the shipped `gate` column).

Same protocol as scripts/calibration_v2.py: band-holdout CV (train 3 ANI
bands, test held-out band, rotate), per-fold median impute + StandardScaler
+ RidgeCV. No flag filtering; only non-finite `ani_gated` rows dropped
(same 21 pairs as v2 -> n = 2,053).

Variants:
  v3a: ani_gated, ani_uniform, af_query, af_reference, std_err, retention,
       n_anchors, n_chains, n_tags           (v2 Set A with ani -> ani_gated)
  v3b: v3a + gate_fallback (0/1 indicator)
Reference: v2 Set A rerun on the gated matrix (ani instead of ani_gated).

External validation: (a) oral/gut 100 same-species pairs vs FastANI, with
ani_gated recomputed in Python from ani_from_loss/ani_from_hist via the
verified gate rule (the oral/gut current matrix predates the gate);
(b) mid-ANI 15 pairs vs ANIm from results/gating_flag/midani_15_gated.tsv
(already gated).
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

import calibration_v2 as c2  # reuse protocol helpers

ROOT = Path(__file__).resolve().parent.parent
GATED = ROOT / "results/anim_truth_2074_gated.tsv"
ACC2SEQ = ROOT / "results/anim_2074_acc2seqid.tsv"
EVAL = ROOT / "results/panel_by_band/eval_pairs.tsv"
ORALGUT_FEAT = ROOT / "results/oral_gut_1225_v8current.tsv"
ORALGUT_MAP = ROOT / "results/oral_gut_1225_acc2seqid.tsv"
ORALGUT_META = ROOT / "data/oral_gut_validation_merged_v8.tsv"
MID_GATED = ROOT / "results/gating_flag/midani_15_gated.tsv"
MID_MAP = ROOT / "results/validation_mid_ani_anim/anim/accession_seqid_map.tsv"
MID_TRUTH = ROOT / "results/validation_mid_ani_anim/anim/anim_truth.tsv"
OUT_CV = ROOT / "results/panel_by_band/calibration_v3_cv.tsv"
OUT_EXT = ROOT / "results/panel_by_band/calibration_v3_external.tsv"
OUT_JSON = ROOT / "results/panel_by_band/linear_cal_v3.json"

V3A = ["ani_gated", "ani_uniform", "af_query", "af_reference", "std_err",
       "retention", "n_anchors", "n_chains", "n_tags"]
V3B = V3A + ["gate_fallback"]
# Rust field order for the LinearCalModel JSON (TSV ani_gated <-
# res.ani_gated*100; ani_uniform <- res.ani*100; std_err <- res.std_err*100;
# gate_fallback <- res.gate_fallback as f64)
RUST_NAMES_V3B = ["ani_gated", "ani_uniform", "af_query", "af_reference",
                  "std_err", "retention", "n_anchors", "n_chains",
                  "n_tags_in_chains", "gate_fallback"]


def load_gated():
    mat = pd.read_csv(GATED, sep="\t")
    acc = pd.read_csv(ACC2SEQ, sep="\t", header=None,
                      names=["accession", "seqid"])
    s2a = dict(zip(acc["seqid"], acc["accession"]))
    mat["query_asm"] = mat["query"].map(s2a)
    mat["ref_asm"] = mat["reference"].map(s2a)
    ep = pd.read_csv(EVAL, sep="\t")
    df = mat.merge(ep[["query_asm", "ref_asm", "band", "anim_ani"]],
                   on=["query_asm", "ref_asm"], how="inner")
    assert len(df) == 2074, len(df)
    df = df.replace([np.inf, -np.inf], np.nan)
    df["gate_fallback"] = (df["gate"] == "uniform_fallback").astype(float)
    n_drop = int(df["ani_gated"].isna().sum())
    df = df[df["ani_gated"].notna()].reset_index(drop=True)
    print(f"gated training rows: {len(df)} (dropped {n_drop} non-finite; "
          f"fallbacks: {int(df['gate_fallback'].sum())})")
    return df


def apply_gate_rule(df):
    """Recompute ani_gated/gate_fallback from loss/hist (verified rule)."""
    df = df.copy()
    fb = (df["ani_from_loss"] - df["ani_from_hist"]).abs() > 5.0
    df["gate_fallback"] = fb.astype(float)
    df["ani_gated"] = np.where(fb, df["ani_uniform"], df["ani"])
    df.loc[df["ani"].isna(), "ani_gated"] = np.nan
    return df


def train_final(df, features, rust_names):
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import RidgeCV
    from sklearn.preprocessing import StandardScaler
    X = df[features].values
    y = df["anim_ani"].values
    imputer = SimpleImputer(strategy="median")
    scaler = StandardScaler()
    Xs = scaler.fit_transform(imputer.fit_transform(X))
    model = RidgeCV(alphas=c2.ALPHAS + [1000], cv=5)
    model.fit(Xs, y)
    mae = float(np.abs(model.predict(Xs) - y).mean())
    cal = {
        "name": "gtdb_r207_linear_cal_v3",
        "feature_names": rust_names,
        "means": scaler.mean_.tolist(),
        "scales": scaler.scale_.tolist(),
        "coefficients": model.coef_.tolist(),
        "intercept": float(model.intercept_),
        "imputer_medians": imputer.statistics_.tolist(),
        "training_n": int(len(df)),
        "training_mae": mae,
    }
    return cal, ("ridge", imputer, scaler, model), mae, float(model.alpha_)


def external_oralgut(fitted, features, tag):
    feat = pd.read_csv(ORALGUT_FEAT, sep="\t")
    amap = pd.read_csv(ORALGUT_MAP, sep="\t", header=None,
                       names=["accession", "seqid"])
    s2a = dict(zip(amap["seqid"], amap["accession"]))
    feat["query_asm"] = feat["query"].map(s2a)
    feat["ref_asm"] = feat["reference"].map(s2a)
    meta = pd.read_csv(ORALGUT_META, sep="\t")
    df = feat.merge(meta[["query", "reference", "label", "q_species",
                          "r_species", "fastani_ani"]],
                    left_on=["query_asm", "ref_asm"],
                    right_on=["query", "reference"], how="inner")
    assert len(df) == 1225, len(df)
    df = apply_gate_rule(df)
    ss = df[(df["q_species"] == df["r_species"]) & df["ani_gated"].notna()
            & df["fastani_ani"].notna()].copy()
    ss["fastani_pct"] = ss["fastani_ani"] * 100.0
    assert len(ss) == 100, len(ss)
    ss["cal"] = c2.predict_dataframe(fitted, ss, features)
    rows = []
    for col, name in [("ani_gated", "raw_gated"), ("cal", tag)]:
        err = ss[col] - ss["fastani_pct"]
        rows.append({"dataset": "oralgut_same_species_vs_fastani_gated",
                     "method": name, "n": len(ss),
                     "MAE": float(err.abs().mean()), "bias": float(err.mean()),
                     "r": float(ss[col].corr(ss["fastani_pct"]))})
    return rows, ss


def external_midani(fitted, features, tag):
    mid = pd.read_csv(MID_GATED, sep="\t")
    amap = pd.read_csv(MID_MAP, sep="\t", header=None,
                       names=["accession", "seqid"])
    s2a = dict(zip(amap["seqid"], amap["accession"]))
    mid["query_asm"] = mid["query"].map(s2a)
    mid["ref_asm"] = mid["reference"].map(s2a)
    mid["gate_fallback"] = (mid["gate"] == "uniform_fallback").astype(float)
    truth = pd.read_csv(MID_TRUTH, sep="\t")
    df = mid.merge(truth, left_on=["query_asm", "ref_asm"],
                   right_on=["query", "reference"], how="inner")
    assert len(df) == 15, len(df)
    df["cal"] = c2.predict_dataframe(fitted, df, features)
    rows = []
    for col, name in [("ani_gated", "raw_gated"), ("cal", tag)]:
        err = df[col] - df["anim_ani"]
        rows.append({"dataset": "midani_15_vs_anim_gated", "method": name,
                     "n": len(df), "MAE": float(err.abs().mean()),
                     "bias": float(err.mean()),
                     "r": float(df[col].corr(df["anim_ani"]))})
    return rows, df


def main():
    df = load_gated()
    y = df["anim_ani"].values
    bands = df["band"].values
    rows = []
    rows += c2.metric_rows("raw_gated", df["ani_gated"].values, y, bands)
    for name, feats in [("v2_setA_rerun", c2.SET_A),
                        ("v3a_gated9", V3A), ("v3b_gated9+gate", V3B)]:
        preds = c2.band_holdout_cv(df, feats, "ridge")
        rows += c2.metric_rows(name, preds, y, bands)
    cv = pd.DataFrame(rows).round(4)
    cv.to_csv(OUT_CV, sep="\t", index=False)
    pd.set_option("display.width", 160)
    print(cv.to_string(index=False))

    def overall(name):
        return cv[(cv.experiment == name) & (cv.band == "all")]["MAE"].iloc[0]

    mae_a, mae_b = overall("v3a_gated9"), overall("v3b_gated9+gate")
    use_b = mae_b < 0.995 * mae_a  # gate indicator must earn its parameter
    feats = V3B if use_b else V3A
    names = RUST_NAMES_V3B if use_b else RUST_NAMES_V3B[:9]
    print(f"\nv3a={mae_a:.4f} v3b={mae_b:.4f} -> final: {'v3b' if use_b else 'v3a'}")

    cal, fitted, mae_in, alpha = train_final(df, feats, names)
    with open(OUT_JSON, "w") as fh:
        json.dump(cal, fh, indent=2)
    print(f"final v3 in-sample MAE {mae_in:.4f}, alpha {alpha}; wrote {OUT_JSON}")

    ext, og = external_oralgut(fitted, feats, f"calibrated_v3")
    mrows, md = external_midani(fitted, feats, f"calibrated_v3")
    ext += mrows
    ext = pd.DataFrame(ext).round(4)
    ext.to_csv(OUT_EXT, sep="\t", index=False)
    print("\n=== external validation (v3) ===")
    print(ext.to_string(index=False))
    print(f"\noral/gut fallbacks recomputed: {int(og['gate_fallback'].sum())}/100; "
          f"mid-ANI fallbacks: {int(md['gate_fallback'].sum())}/15")


if __name__ == "__main__":
    main()
