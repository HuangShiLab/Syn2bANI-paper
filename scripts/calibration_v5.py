#!/usr/bin/env python3
"""Calibration v5: v4 protocol retrained on POST-RESCUE feature matrices.

The short-contig rescue pass (Syn2bANI 0aabd0c) changed estimator features
(af_query, retention, n_tags_in_chains, ...) on fragmented inputs, so the
deployed v4 model (trained on pre-rescue features) is inconsistent with the
shipped estimator on draft genomes. v5 retrains on matrices re-run with the
068119c binary:
  results/anim_truth_2074_v9rescue.tsv   (2,074 GTDB pairs)
  results/anim_truth_hi95_v9rescue.tsv   (467 hi95 pairs)
Externals are also post-rescue re-runs:
  results/oral_gut_1225_v9rescue.tsv          (vs FastANI, 100 same-species)
  results/gating_flag/midani_15_v9rescue.tsv  (vs ANIm, 15 pairs)

Protocol identical to v4 (scripts/calibration_v4.py): band-holdout CV
(train 3 bands, test held-out, rotate), per-fold median impute +
StandardScaler + RidgeCV, training rows shuffled before CV, seeds 0-4 with
seed 0 primary.

Variants:
  v5a: V3A features (ani_gated, ani_uniform, af_query, af_reference,
       std_err, retention, n_anchors, n_chains, n_tags)
  v5b: v5a + gate_fallback (0/1 indicator)
Selection rule (same as v3/v4): v5b only if overall CV MAE < 0.995 * v5a.
v4a/v4b are rerun on the old 2,053-row subset of the POST-RESCUE 2074
matrix with the same seeded protocol as the reference line.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

import calibration_v2 as c2  # protocol helpers (band_holdout_cv, metric_rows)
import calibration_v3 as c3  # feature sets / rust names / gate rule
import calibration_v4 as c4  # shuffled(), SEEDS, PRIMARY_SEED, band constants

ROOT = Path(__file__).resolve().parent.parent
GATED = ROOT / "results/anim_truth_2074_v9rescue.tsv"
ACC2SEQ = ROOT / "results/anim_2074_acc2seqid.tsv"
EVAL = ROOT / "results/panel_by_band/eval_pairs.tsv"
HI95_GATED = ROOT / "results/anim_truth_hi95_v9rescue.tsv"
HI95_TRUTH = ROOT / "results/anim_truth_hi95.tsv"
ORALGUT_FEAT = ROOT / "results/oral_gut_1225_v9rescue.tsv"
ORALGUT_MAP = ROOT / "results/oral_gut_1225_acc2seqid.tsv"
ORALGUT_META = ROOT / "data/oral_gut_validation_merged_v8.tsv"
MID_GATED = ROOT / "results/gating_flag/midani_15_v9rescue.tsv"
MID_MAP = ROOT / "results/validation_mid_ani_anim/anim/accession_seqid_map.tsv"
MID_TRUTH = ROOT / "results/validation_mid_ani_anim/anim/anim_truth.tsv"
OUT_CV = ROOT / "results/panel_by_band/calibration_v5_cv.tsv"
OUT_EXT = ROOT / "results/panel_by_band/calibration_v5_external.tsv"
OUT_JSON = ROOT / "results/panel_by_band/linear_cal_v5.json"

SEEDS = c4.SEEDS
PRIMARY_SEED = c4.PRIMARY_SEED
TOP_BAND = c4.TOP_BAND


def load_gated_v5():
    """Post-rescue 2074 matrix; identical join/asserts as c3.load_gated."""
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
    print(f"gated v9rescue training rows: {len(df)} (dropped {n_drop} "
          f"non-finite; fallbacks: {int(df['gate_fallback'].sum())})")
    return df


def load_hi95_v5():
    """Post-rescue hi95 matrix; identical to c4.load_hi95."""
    mat = pd.read_csv(HI95_GATED, sep="\t")
    truth = pd.read_csv(HI95_TRUTH, sep="\t")
    df = mat.merge(truth[["query", "reference", "anim_ani"]],
                   on=["query", "reference"], how="inner")
    assert len(df) == 467, len(df)
    in_band = df["anim_ani"].between(c4.HI95_LO, c4.HI95_HI, inclusive="left")
    assert in_band.all(), "hi95 pairs outside the 95-99.5 band"
    df["band"] = TOP_BAND
    df = df.replace([np.inf, -np.inf], np.nan)
    df["gate_fallback"] = (df["gate"] == "uniform_fallback").astype(float)
    n_drop = int(df["ani_gated"].isna().sum())
    df = df[df["ani_gated"].notna()].reset_index(drop=True)
    print(f"hi95 v9rescue training rows: {len(df)} (dropped {n_drop} "
          f"non-finite; fallbacks: {int(df['gate_fallback'].sum())})")
    return df


def verify_gate_rule(df, tag):
    """The shipped gate column must match the verified Python rule."""
    chk = c3.apply_gate_rule(df)
    same = (chk["gate_fallback"].values == df["gate_fallback"].values)
    gated_eq = np.allclose(chk["ani_gated"].values, df["ani_gated"].values,
                           equal_nan=True)
    print(f"{tag}: shipped gate matches Python rule on "
          f"{int(same.sum())}/{len(df)} rows; ani_gated equal: {gated_eq}")
    assert same.all() and gated_eq, f"gate rule mismatch on {tag}"


def external_oralgut_v5(fitted, features, tag):
    """c3.external_oralgut on the post-rescue 1225-pair matrix."""
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
    df["gate_fallback"] = (df["gate"] == "uniform_fallback").astype(float)
    verify_gate_rule(df, "oral/gut v9rescue")
    df = c3.apply_gate_rule(df)  # identical to shipped gate (verified above)
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


def external_midani_v5(fitted, features, tag):
    """c3.external_midani on the post-rescue mid-ANI matrix."""
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


def train_final_v5(df, features, rust_names, seed):
    cal, fitted, mae, alpha = c4.train_final(df, features, rust_names, seed)
    cal["name"] = "gtdb_r207_linear_cal_v5"
    return cal, fitted, mae, alpha


def main():
    old = load_gated_v5()        # post-rescue 2074 subset, bands from eval_pairs
    new = load_hi95_v5()         # post-rescue 467 rows, band = 0.95-0.99
    verify_gate_rule(old, "2074 v9rescue")
    verify_gate_rule(new, "hi95 v9rescue")
    df = pd.concat([old, new], ignore_index=True)
    print(f"merged training rows: {len(df)} "
          f"(top band n: {int((old.band == TOP_BAND).sum())} -> "
          f"{int((df.band == TOP_BAND).sum())})")

    rows = []
    # raw reference on the merged set
    d0 = c4.shuffled(df, PRIMARY_SEED)
    for r in c2.metric_rows("raw_gated", d0["ani_gated"].values,
                            d0["anim_ani"].values, d0["band"].values):
        r["seed"] = PRIMARY_SEED
        rows.append(r)
    # v4 reference rerun on the old 2,053-row subset, same seeded protocol
    rows += c4.cv_experiment(old, "v4a_ref_gated9", c3.V3A, SEEDS)
    rows += c4.cv_experiment(old, "v4b_ref_gated9+gate", c3.V3B, SEEDS)
    # v5 variants on the merged 2,520-row set
    rows += c4.cv_experiment(df, "v5a_gated9", c3.V3A, SEEDS)
    rows += c4.cv_experiment(df, "v5b_gated9+gate", c3.V3B, SEEDS)

    cv = pd.DataFrame(rows)
    cv = c4.add_mean_rows(cv).round(4)
    cv = cv[["experiment", "seed", "band", "n", "MAE", "RMSE", "bias", "r"]]
    cv.to_csv(OUT_CV, sep="\t", index=False)
    pd.set_option("display.width", 160)
    print(cv[cv.band == "all"].to_string(index=False))

    print("\n=== seed sensitivity (overall CV MAE) ===")
    for name in ["v4a_ref_gated9", "v4b_ref_gated9+gate", "v5a_gated9",
                 "v5b_gated9+gate"]:
        per_seed = [c4.overall(cv, name, s) for s in SEEDS]
        print(f"{name:22s} " +
              " ".join(f"{v:.4f}" for v in per_seed) +
              f"  mean={np.mean(per_seed):.4f} spread="
              f"{max(per_seed) - min(per_seed):.4f}")

    mae_a = c4.overall(cv, "v5a_gated9")
    mae_b = c4.overall(cv, "v5b_gated9+gate")
    use_b = mae_b < 0.995 * mae_a  # gate indicator must earn its parameter
    decisions = {s: c4.overall(cv, "v5b_gated9+gate", s)
                 < 0.995 * c4.overall(cv, "v5a_gated9", s) for s in SEEDS}
    feats = c3.V3B if use_b else c3.V3A
    names = c3.RUST_NAMES_V3B if use_b else c3.RUST_NAMES_V3B[:9]
    tag = "v5b" if use_b else "v5a"
    print(f"\nv5a={mae_a:.4f} v5b={mae_b:.4f} (seed {PRIMARY_SEED}) "
          f"-> final: {tag}; use_b per seed: {decisions}")

    cal, fitted, mae_in, alpha = train_final_v5(df, feats, names, PRIMARY_SEED)
    with open(OUT_JSON, "w") as fh:
        json.dump(cal, fh, indent=2)
    print(f"final v5 in-sample MAE {mae_in:.4f}, alpha {alpha}; "
          f"wrote {OUT_JSON}")

    ext, og = external_oralgut_v5(fitted, feats, "calibrated_v5")
    mrows, md = external_midani_v5(fitted, feats, "calibrated_v5")
    ext += mrows
    ext = pd.DataFrame(ext).round(4)
    ext.to_csv(OUT_EXT, sep="\t", index=False)
    print("\n=== external validation (v5) ===")
    print(ext.to_string(index=False))
    print(f"\noral/gut fallbacks recomputed: "
          f"{int(og['gate_fallback'].sum())}/100; mid-ANI fallbacks: "
          f"{int(md['gate_fallback'].sum())}/15")

    # gates vs deployed v4 (primary seed, per-band v5-selected vs v4a_ref
    # rerun on the same post-rescue old subset; external gates vs the v4
    # external numbers, which were computed on PRE-rescue features)
    print("\n=== gates vs v4 ===")
    v4_ref = "v4a_ref_gated9"
    v5_sel = "v5b_gated9+gate" if use_b else "v5a_gated9"
    ok = True
    for band in c2.BAND_ORDER:
        if band == TOP_BAND:
            continue
        m4 = cv[(cv.experiment == v4_ref) & (cv.band == band) &
                (cv.seed == PRIMARY_SEED)]["MAE"].iloc[0]
        m5 = cv[(cv.experiment == v5_sel) & (cv.band == band) &
                (cv.seed == PRIMARY_SEED)]["MAE"].iloc[0]
        passed = m5 <= m4 + 0.02
        ok &= passed
        print(f"band {band}: v4_ref {m4:.4f} -> v5 {m5:.4f} "
              f"(delta {m5 - m4:+.4f}) {'PASS' if passed else 'FAIL'}")
    og5 = ext[(ext.dataset == "oralgut_same_species_vs_fastani_gated") &
              (ext.method == "calibrated_v5")]["MAE"].iloc[0]
    mid5 = ext[(ext.dataset == "midani_15_vs_anim_gated") &
               (ext.method == "calibrated_v5")]["MAE"].iloc[0]
    og_ok = og5 <= 0.4634   # v4 external (pre-rescue features)
    mid_ok = mid5 <= 0.7942
    print(f"oral/gut external: v4 0.4634 -> v5 {og5:.4f} "
          f"{'PASS' if og_ok else 'FAIL'}")
    print(f"mid-ANI external:  v4 0.7942 -> v5 {mid5:.4f} "
          f"{'PASS' if mid_ok else 'FAIL'}")
    ok &= og_ok and mid_ok
    print(f"ALL GATES {'PASS' if ok else 'FAIL'}")


if __name__ == "__main__":
    main()
