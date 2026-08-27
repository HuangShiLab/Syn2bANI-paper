#!/usr/bin/env python3
"""Calibration v2: retrain on current-binary (v8 @ 69ce9f4) features and
validate externally.

Protocol identical to scripts/anim_main_table_4e.py: band-holdout CV (train
on 3 ANI bands, test on held-out band, rotate), per-fold median imputation +
StandardScaler + RidgeCV (alphas 1e-3..100, inner cv=5). No flag filtering
(the ok/INCONSISTENT flag is inverted on GTDB); only rows with non-finite
`ani` are dropped (21 of 2,074). inf feature values (het_shape) are treated
as missing and median-imputed, matching the Rust predict() non-finite ->
imputer-median path.

Feature sets:
  A (base 9):  ani, ani_uniform, af_query, af_reference, std_err, retention,
               n_anchors, n_chains, n_tags
  B (expanded): A + anchor_adjacency, breakpoint_count, enzyme_spread,
               enzyme_chi2, het_shape, ani_from_loss, ani_from_hist,
               max_block_anchors, mean_block_anchors
  C: set B with GradientBoostingRegressor (sklearn defaults) instead of ridge.

External validation: (a) oral/gut same-species 100 pairs vs FastANI
(Set-A model only — the file lacks the synteny/dispersion columns);
(b) 15 mid-ANI pairs vs ANIm truth (all Set-B columns present).
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent.parent
MATRIX = ROOT / "results/anim_truth_2074_v8current.tsv"
ACC2SEQ = ROOT / "results/anim_2074_acc2seqid.tsv"
EVAL = ROOT / "results/panel_by_band/eval_pairs.tsv"
ORALGUT = ROOT / "data/oral_gut_validation_merged_v8.tsv"
MID_S2B = ROOT / "results/validation_mid_ani_anim/anim/syn2bani_v8_4e.tsv"
MID_MAP = ROOT / "results/validation_mid_ani_anim/anim/accession_seqid_map.tsv"
MID_TRUTH = ROOT / "results/validation_mid_ani_anim/anim/anim_truth.tsv"
OUT_CV = ROOT / "results/panel_by_band/calibration_v2_cv.tsv"
OUT_EXT = ROOT / "results/panel_by_band/calibration_v2_external.tsv"
OUT_JSON = ROOT / "results/panel_by_band/linear_cal_v2.json"

ALPHAS = [1e-3, 1e-2, 0.1, 1, 10, 100]
BAND_ORDER = ["0.8-0.85", "0.85-0.9", "0.9-0.95", "0.95-0.99"]

SET_A = ["ani", "ani_uniform", "af_query", "af_reference", "std_err",
         "retention", "n_anchors", "n_chains", "n_tags"]
SET_B_EXTRA = ["anchor_adjacency", "breakpoint_count", "enzyme_spread",
               "enzyme_chi2", "het_shape", "ani_from_loss", "ani_from_hist",
               "max_block_anchors", "mean_block_anchors"]
SET_B = SET_A + SET_B_EXTRA

# LinearCalModel feature names, in the order Rust's predict_from_result
# supplies them (TSV ani <- res.ani_het*100; TSV ani_uniform <- res.ani*100;
# TSV std_err <- res.std_err*100; TSV enzyme_spread <- agreement.spread*100;
# TSV ani_from_loss/hist <- fields *100; all others direct).
RUST_NAMES = ["ani_het", "ani_uniform", "af_query", "af_reference", "std_err",
              "retention", "n_anchors", "n_chains", "n_tags_in_chains",
              "anchor_adjacency", "breakpoint_count", "enzyme_spread",
              "enzyme_chi2", "het_shape", "ani_from_loss", "ani_from_hist",
              "max_block_anchors", "mean_block_anchors"]


def load_training():
    mat = pd.read_csv(MATRIX, sep="\t")
    acc = pd.read_csv(ACC2SEQ, sep="\t", header=None,
                      names=["accession", "seqid"])
    s2a = dict(zip(acc["seqid"], acc["accession"]))
    mat["query_asm"] = mat["query"].map(s2a)
    mat["ref_asm"] = mat["reference"].map(s2a)
    ep = pd.read_csv(EVAL, sep="\t")
    df = mat.merge(ep[["query_asm", "ref_asm", "band", "anim_ani"]],
                   on=["query_asm", "ref_asm"], how="inner")
    assert len(df) == 2074, len(df)
    # non-finite -> NaN so the imputer sees them as missing (het_shape=inf)
    df = df.replace([np.inf, -np.inf], np.nan)
    n_drop = int(df["ani"].isna().sum())
    df = df[df["ani"].notna()].reset_index(drop=True)
    print(f"training rows: {len(df)} (dropped {n_drop} non-finite ani; "
          f"no flag filtering)")
    return df


def make_model(kind, Xtr, ytr):
    imputer = SimpleImputer(strategy="median")
    Xi = imputer.fit_transform(Xtr)
    if kind == "ridge":
        scaler = StandardScaler()
        Xs = scaler.fit_transform(Xi)
        model = RidgeCV(alphas=ALPHAS, cv=5)
        model.fit(Xs, ytr)
        return ("ridge", imputer, scaler, model)
    model = GradientBoostingRegressor(random_state=0)
    model.fit(Xi, ytr)
    return ("gbrt", imputer, None, model)


def model_predict(fitted, X):
    kind, imputer, scaler, model = fitted
    Xi = imputer.transform(X)
    if scaler is not None:
        Xi = scaler.transform(Xi)
    return model.predict(Xi)


def band_holdout_cv(df, features, kind):
    X = df[features].values
    y = df["anim_ani"].values
    bands = df["band"].values
    preds = np.full(len(df), np.nan)
    for band in BAND_ORDER:
        test = bands == band
        fitted = make_model(kind, X[~test], y[~test])
        preds[test] = model_predict(fitted, X[test])
    return preds


def metric_rows(name, pred, truth, bands):
    rows = []
    for band in BAND_ORDER + ["all"]:
        m = np.isfinite(pred) if band == "all" else (bands == band) & np.isfinite(pred)
        if m.sum() < 5:
            continue
        err = pred[m] - truth[m]
        rows.append({"experiment": name, "band": band, "n": int(m.sum()),
                     "MAE": float(np.abs(err).mean()),
                     "RMSE": float(np.sqrt((err ** 2).mean())),
                     "bias": float(err.mean()),
                     "r": float(np.corrcoef(pred[m], truth[m])[0, 1])})
    return rows


def train_final_ridge(df, features):
    X = df[features].values
    y = df["anim_ani"].values
    imputer = SimpleImputer(strategy="median")
    scaler = StandardScaler()
    Xs = scaler.fit_transform(imputer.fit_transform(X))
    model = RidgeCV(alphas=ALPHAS + [1000], cv=5)
    model.fit(Xs, y)
    mae = float(np.abs(model.predict(Xs) - y).mean())
    names = [RUST_NAMES[SET_B.index(f)] if f in SET_B else f for f in features]
    cal = {
        "name": "gtdb_r207_linear_cal_v2",
        "feature_names": names,
        "means": scaler.mean_.tolist(),
        "scales": scaler.scale_.tolist(),
        "coefficients": model.coef_.tolist(),
        "intercept": float(model.intercept_),
        "imputer_medians": imputer.statistics_.tolist(),
        "training_n": int(len(df)),
        "training_mae": mae,
    }
    return cal, ("ridge", imputer, scaler, model), mae, float(model.alpha_)


def predict_dataframe(fitted, df, features):
    X = df[features].replace([np.inf, -np.inf], np.nan).values
    return model_predict(fitted, X)


def external_oralgut(fitted_a):
    og = pd.read_csv(ORALGUT, sep="\t")
    ss = og[og["q_species"] == og["r_species"]].copy()
    ss["fastani_pct"] = ss["fastani_ani"] * 100.0
    ss = ss[ss["ani"].notna() & ss["fastani_pct"].notna()]
    ss["cal"] = predict_dataframe(fitted_a, ss, SET_A)
    rows = []
    for col, name in [("ani", "raw_gamma"), ("ani_uniform", "raw_uniform"),
                      ("cal", "calibrated_v2_setA")]:
        err = ss[col] - ss["fastani_pct"]
        rows.append({"dataset": "oralgut_same_species_vs_fastani",
                     "method": name, "n": len(ss),
                     "MAE": float(err.abs().mean()),
                     "bias": float(err.mean()),
                     "r": float(ss[col].corr(ss["fastani_pct"]))})
    return rows


def external_midani(fitted, features, tag):
    mid = pd.read_csv(MID_S2B, sep="\t")
    amap = pd.read_csv(MID_MAP, sep="\t", header=None,
                       names=["accession", "seqid"])
    s2a = dict(zip(amap["seqid"], amap["accession"]))
    mid["query_asm"] = mid["query"].map(s2a)
    mid["ref_asm"] = mid["reference"].map(s2a)
    truth = pd.read_csv(MID_TRUTH, sep="\t")
    df = mid.merge(truth, left_on=["query_asm", "ref_asm"],
                   right_on=["query", "reference"], how="inner")
    assert len(df) == 15, len(df)
    df["cal"] = predict_dataframe(fitted, df, features)
    rows = []
    for col, name in [("ani", "raw_gamma"), ("ani_uniform", "raw_uniform"),
                      ("cal", f"calibrated_v2_{tag}")]:
        err = df[col] - df["anim_ani"]
        rows.append({"dataset": "midani_15_vs_anim", "method": name,
                     "n": len(df), "MAE": float(err.abs().mean()),
                     "bias": float(err.mean()),
                     "r": float(df[col].corr(df["anim_ani"]))})
    return rows, df[["query_asm", "ref_asm", "anim_ani", "ani",
                     "ani_uniform", "cal", "flag"]]


def main():
    df = load_training()
    y = df["anim_ani"].values
    bands = df["band"].values
    rows = []
    # raw references
    rows += metric_rows("raw_gamma", df["ani"].values, y, bands)
    for name, feats, kind in [("A_base9_ridge", SET_A, "ridge"),
                              ("B_expanded_ridge", SET_B, "ridge"),
                              ("C_expanded_gbrt", SET_B, "gbrt")]:
        preds = band_holdout_cv(df, feats, kind)
        rows += metric_rows(name, preds, y, bands)
    cv = pd.DataFrame(rows).round(4)
    cv.to_csv(OUT_CV, sep="\t", index=False)
    pd.set_option("display.width", 160)
    print(cv.to_string(index=False))

    mae_a = cv[(cv.experiment == "A_base9_ridge") & (cv.band == "all")]["MAE"].iloc[0]
    mae_b = cv[(cv.experiment == "B_expanded_ridge") & (cv.band == "all")]["MAE"].iloc[0]
    # B must beat A by >2% relative to count as a clear win
    use_b = mae_b < 0.98 * mae_a
    final_feats = SET_B if use_b else SET_A
    print(f"\nCV all-band MAE: A={mae_a:.4f} B={mae_b:.4f} "
          f"-> final feature set: {'B' if use_b else 'A'}")

    cal_a, fitted_a, mae_a_in, alpha_a = train_final_ridge(df, SET_A)
    cal_b, fitted_b, mae_b_in, alpha_b = train_final_ridge(df, SET_B)
    final_cal = cal_b if use_b else cal_a
    final_fitted = fitted_b if use_b else fitted_a
    with open(OUT_JSON, "w") as fh:
        json.dump(final_cal, fh, indent=2)
    print(f"final model: {'B' if use_b else 'A'}, in-sample MAE "
          f"{mae_b_in if use_b else mae_a_in:.4f}, alpha "
          f"{alpha_b if use_b else alpha_a}, wrote {OUT_JSON}")
    # also keep the Set-A model JSON for the oral/gut fallback path
    with open(ROOT / "results/panel_by_band/linear_cal_v2_setA.json", "w") as fh:
        json.dump(cal_a, fh, indent=2)

    ext = external_oralgut(fitted_a)
    mid_rows, mid_detail = external_midani(final_fitted, final_feats,
                                           "setB" if use_b else "setA")
    ext += mid_rows
    ext = pd.DataFrame(ext).round(4)
    ext.to_csv(OUT_EXT, sep="\t", index=False)
    mid_detail.round(4).to_csv(
        ROOT / "results/panel_by_band/calibration_v2_midani_pairs.tsv",
        sep="\t", index=False)
    print("\n=== external validation ===")
    print(ext.to_string(index=False))
    print("\n=== mid-ANI per-pair ===")
    print(mid_detail.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
