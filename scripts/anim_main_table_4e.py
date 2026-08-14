#!/usr/bin/env python3
"""Consolidated ANIm main table + band-holdout ridge calibration on the
current 4-enzyme (BcgI,AlfI,AloI,FalI) feature columns.

Methodology follows scripts/analyze_error_drivers.py (calibration_experiment):
median imputation + RidgeCV, band-holdout CV (train on 3 bands, test on the
held-out band, rotate). Difference vs the original: a StandardScaler is fitted
on each training fold and applied to the held-out band (the original CV probe
did not standardize; only its final serialized model did). Imputer and scaler
are both fitted on the training folds only.

BELOW_DETECTION handling: pairs whose current 4e s2b_ani is NaN
(flag == BELOW_DETECTION) are EXCLUDED from the 4e feature-based methods
(gamma, uniform, ridge_cv) and n is reported; the ridge needs s2b_ani as its
primary feature and these pairs carry no point estimate. They remain in the
11e / skani rows (those estimators produced values for all 2074 pairs).
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parent.parent
EVAL = ROOT / "results/panel_by_band/eval_pairs.tsv"
MATRIX = ROOT / "results/matrix_gtdb_r207_100k_v8_final.tsv"
OUT_TABLE = ROOT / "results/panel_by_band/anim_main_table.tsv"
OUT_PREDS = ROOT / "results/panel_by_band/ridge_cv_preds_4e.tsv"
OUT_JSON = ROOT / "results/panel_by_band/ridge_cv_4e_report.json"

FEATURES = ["s2b_ani", "s2b_ani_uniform", "s2b_af_q", "s2b_af_r",
            "s2b_std_err", "s2b_retention", "s2b_n_anchors",
            "s2b_n_chains", "s2b_n_tags"]
ALPHAS = [1e-3, 1e-2, 0.1, 1, 10, 100]
BAND_ORDER = ["0.8-0.85", "0.85-0.9", "0.9-0.95", "0.95-0.99"]


def load():
    ep = pd.read_csv(EVAL, sep="\t")
    mat = pd.read_csv(MATRIX, sep="\t")
    mat = mat.rename(columns={"query": "query_asm", "reference": "ref_asm"})
    df = ep.merge(mat, on=["query_asm", "ref_asm"], how="left",
                  suffixes=("", "_mat"))
    assert len(df) == len(ep), (len(df), len(ep))
    # matrix skani/fastani are fractions 0-1 -> percent
    df["fastani_pct"] = df["fastani_ani"] * 100.0
    return df


def ridge_band_holdout(df):
    """Return per-pair out-of-fold ridge predictions (NaN where s2b_ani NaN)."""
    usable = df[df["s2b_ani_mat"].notna()].copy()
    # after merge, only s2b_ani collides with the eval_pairs column -> _mat
    Xdf = usable[[c + "_mat" if c == "s2b_ani" else c for c in FEATURES]]
    Xdf.columns = FEATURES
    y = usable["anim_ani"].values
    bands = usable["band"].values
    preds = np.full(len(usable), np.nan)
    coef_sum = np.zeros(len(FEATURES))
    n_fits = 0
    fold_alphas = {}
    for band in BAND_ORDER:
        test = bands == band
        train = ~test
        imputer = SimpleImputer(strategy="median")
        scaler = StandardScaler()
        Xtr = scaler.fit_transform(imputer.fit_transform(Xdf.values[train]))
        Xte = scaler.transform(imputer.transform(Xdf.values[test]))
        model = RidgeCV(alphas=ALPHAS, cv=5)
        model.fit(Xtr, y[train])
        preds[test] = model.predict(Xte)
        coef_sum += model.coef_
        n_fits += 1
        fold_alphas[band] = float(model.alpha_)
    usable = usable.copy()
    usable["ridge_pred"] = preds
    return usable, coef_sum / n_fits, fold_alphas


def metrics(err, pred, truth):
    err = np.asarray(err, dtype=float)
    return {
        "n": int(np.isfinite(err).sum()),
        "MAE": float(np.nanmean(np.abs(err))),
        "RMSE": float(np.sqrt(np.nanmean(err ** 2))),
        "bias": float(np.nanmean(err)),
        "r": float(np.corrcoef(pred, truth)[0, 1]),
    }


def build_table(df, ridge):
    rows = []

    def add(method, pred, truth, bands, subset_label=None):
        pred = np.asarray(pred, dtype=float)
        truth = np.asarray(truth, dtype=float)
        valid = np.isfinite(pred) & np.isfinite(truth)
        for band in list(BAND_ORDER) + ["all"]:
            m = valid if band == "all" else valid & (bands == band)
            if m.sum() < 5:
                continue
            met = metrics(pred[m] - truth[m], pred[m], truth[m])
            rows.append({"method": method, "band": band, **met})

    bands = df["band"].values
    truth = df["anim_ani"].values
    add("syn2bani_4e_gamma", df["s2b_ani_mat"].values, truth, bands)
    add("syn2bani_4e_uniform", df["s2b_ani_uniform"].values, truth, bands)
    rb = ridge["band"].values
    add("syn2bani_4e_ridge_cv", ridge["ridge_pred"].values,
        ridge["anim_ani"].values, rb)
    add("syn2bani_11e", df["s2b_ani"].values, truth, bands)
    add("syn2bani_11e_cal", df["s2b_ani_cal"].values, truth, bands)
    add("skani", df["skani_ani"].values, truth, bands)
    fa = df["fastani_pct"].values
    add("FastANI_subset", fa, truth, bands)
    return pd.DataFrame(rows)


def main():
    df = load()
    n_below = int(df["s2b_ani_mat"].isna().sum())
    flags = df.loc[df["s2b_ani_mat"].isna(), "s2b_flag"].value_counts().to_dict()

    ridge, mean_coef, fold_alphas = ridge_band_holdout(df)
    ridge_out = ridge[["query_asm", "ref_asm", "band", "anim_ani",
                       "s2b_ani_mat", "ridge_pred"]].rename(
        columns={"s2b_ani_mat": "s2b_ani_4e"})
    ridge_out.to_csv(OUT_PREDS, sep="\t", index=False)

    table = build_table(df, ridge)
    table = table.round(4)
    table.to_csv(OUT_TABLE, sep="\t", index=False)

    report = {
        "n_pairs": int(len(df)),
        "n_below_detection_excluded_4e": n_below,
        "below_detection_flags": flags,
        "features": FEATURES,
        "alphas_grid": ALPHAS,
        "fold_alphas": fold_alphas,
        "mean_coef": dict(zip(FEATURES, mean_coef.tolist())),
        "ridge_overall": table[(table.method == "syn2bani_4e_ridge_cv")
                               & (table.band == "all")].to_dict("records"),
        "fastani_subset_n": int(df["fastani_pct"].notna().sum()),
    }
    with open(OUT_JSON, "w") as fh:
        json.dump(report, fh, indent=2)

    pd.set_option("display.width", 160)
    print(table.to_string(index=False))
    print(f"\nexcluded BELOW_DETECTION (4e NaN): {n_below}, flags={flags}")
    print(f"fold alphas: {fold_alphas}")
    print(f"mean coef: {report['mean_coef']}")
    print(f"wrote {OUT_TABLE}\nwrote {OUT_PREDS}\nwrote {OUT_JSON}")


if __name__ == "__main__":
    main()
