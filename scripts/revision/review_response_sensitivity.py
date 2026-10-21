#!/usr/bin/env python3
"""Generate review-response sensitivity tables from archived result tables."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/revision"


def metrics(y: pd.Series, pred: pd.Series) -> tuple[int, float, float, float]:
    ok = y.notna() & pred.notna()
    yy = y[ok].astype(float)
    pp = pred[ok].astype(float)
    if len(yy) == 0:
        return 0, np.nan, np.nan, np.nan
    error = pp - yy
    r = float(np.corrcoef(pp, yy)[0, 1]) if np.std(pp) > 0 and np.std(yy) > 0 else np.nan
    return int(len(yy)), float(np.abs(error).mean()), float(error.mean()), r


def af_benchmark() -> None:
    truth = pd.read_csv(ROOT / "results/gtdb50k/truth_50k.tsv", sep="\t")
    pred = pd.read_csv(ROOT / "results/gtdb50k/s2b_50k.tsv", sep="\t")
    pairs = pd.read_csv(ROOT / "results/gtdb50k/pairs_50k.tsv", sep="\t")
    fast = pd.read_csv(ROOT / "results/gtdb50k/fastani_50k.tsv", sep="\t")
    pairs["pairid"] = pairs["q_acc"] + "__" + pairs["r_acc"]
    df = (
        truth.merge(pred[["pairid", "ani_gated", "ani_cal"]], on="pairid", how="inner")
        .merge(pairs[["pairid", "skani_ani"]], on="pairid", how="left")
        .merge(fast[["pairid", "fastani_ani"]], on="pairid", how="left")
    )
    methods = {
        "syn2bani_raw_gated": "ani_gated",
        "syn2bani_calibrated_v6": "ani_cal",
        "skani": "skani_ani",
        "fastani": "fastani_ani",
    }
    rows = []
    tiers = ["strict", "low-AF", "verylow-AF"]
    for tier in tiers + ["all"]:
        sub = df if tier == "all" else df[df["af_tier"] == tier]
        complete = sub[
            sub[["ani_gated", "ani_cal", "skani_ani", "fastani_ani"]].notna().all(axis=1)
        ]
        for method, col in methods.items():
            n, mae, bias, r = metrics(complete["anim_ani"], complete[col])
            rows.append({
                "af_tier": tier,
                "denominator": "common_complete_case",
                "method": method,
                "n_pairs": n,
                "mae_ani_points": mae,
                "bias_ani_points": bias,
                "pearson_r": r,
            })
        rows.append({
            "af_tier": tier,
            "denominator": "all_pairs_in_stratum",
            "method": "all_methods_complete",
            "n_pairs": len(complete),
            "mae_ani_points": np.nan,
            "bias_ani_points": np.nan,
            "pearson_r": np.nan,
        })
    pd.DataFrame(rows).to_csv(OUT / "af_stratified_benchmark.tsv", sep="\t", index=False, float_format="%.6g")


def gtdb_scope() -> None:
    joint = pd.read_csv(ROOT / "results/census/gtdb/gtdb_species_ani_sv_joint.tsv", sep="\t")
    joint = joint[joint["expected_pairs"].notna()].copy()
    joint["coverage"] = joint["structural_rows"] / joint["expected_pairs"]
    joint["fraction_ge2_reported"] = joint["sv_ge2_pairs"] / joint["structural_rows"]
    taxa = joint["taxonomy"].str.split(";", expand=True)
    joint["genus"] = taxa[4].str.strip() if taxa.shape[1] > 4 else np.nan
    joint["family"] = taxa[3].str.strip() if taxa.shape[1] > 3 else np.nan

    def scope_mask(name: str) -> pd.Series:
        base = joint["expected_pairs"] >= 100
        if name == "complete_joint_expected_ge100":
            return base & (joint["ani_sv_layer"] == "complete")
        if name == "all_joint_expected_ge100":
            return base
        if name == "all_joint_coverage_ge0.99_expected_ge100":
            return base & (joint["coverage"] >= 0.99)
        if name == "all_joint_coverage_ge0.95_expected_ge100":
            return base & (joint["coverage"] >= 0.95)
        if name == "complete_joint_exclude_Escherichia_coli":
            return base & (joint["ani_sv_layer"] == "complete") & (~joint["species_cluster"].str.contains("Escherichia coli", regex=False))
        raise KeyError(name)

    def band(x: float) -> str:
        if pd.isna(x):
            return "missing"
        if x < 95:
            return "<95"
        if x < 97:
            return "95-97"
        if x < 98:
            return "97-98"
        if x < 99:
            return "98-99"
        return ">=99"

    joint["ani_band"] = joint["mean_skani_ani"].map(band)
    scopes = [
        "complete_joint_expected_ge100",
        "all_joint_expected_ge100",
        "all_joint_coverage_ge0.99_expected_ge100",
        "all_joint_coverage_ge0.95_expected_ge100",
        "complete_joint_exclude_Escherichia_coli",
    ]
    rows = []
    for scope in scopes:
        mask = scope_mask(scope)
        sub = joint[mask]
        genus = sub.groupby("genus").apply(
            lambda g: pd.Series({
                "n_species": len(g),
                "pairs": g["structural_rows"].sum(),
                "weighted_fraction": g["sv_ge2_pairs"].sum() / g["structural_rows"].sum(),
            }),
            include_groups=False,
        )
        family = sub.groupby("family").apply(
            lambda g: pd.Series({
                "n_species": len(g),
                "pairs": g["structural_rows"].sum(),
                "weighted_fraction": g["sv_ge2_pairs"].sum() / g["structural_rows"].sum(),
            }),
            include_groups=False,
        )
        rows.append({
            "scope": scope,
            "summary": "all",
            "ani_band": "all",
            "n_species": len(sub),
            "n_genera": int(genus.shape[0]),
            "n_families": int(family.shape[0]),
            "pairs": int(sub["structural_rows"].sum()),
            "unweighted_median_fraction_ge2": float(sub["fraction_ge2_reported"].median()),
            "pair_weighted_fraction_ge2": float(sub["sv_ge2_pairs"].sum() / sub["structural_rows"].sum()),
            "median_genus_weighted_fraction": float(genus["weighted_fraction"].median()),
            "median_family_weighted_fraction": float(family["weighted_fraction"].median()),
        })
        for band_name in ["<95", "95-97", "97-98", "98-99", ">=99"]:
            bsub = sub[sub["ani_band"] == band_name]
            if len(bsub) == 0:
                continue
            bg = bsub.groupby("genus").apply(
                lambda g: pd.Series({
                    "pairs": g["structural_rows"].sum(),
                    "weighted_fraction": g["sv_ge2_pairs"].sum() / g["structural_rows"].sum(),
                }),
                include_groups=False,
            )
            bf = bsub.groupby("family").apply(
                lambda g: pd.Series({
                    "pairs": g["structural_rows"].sum(),
                    "weighted_fraction": g["sv_ge2_pairs"].sum() / g["structural_rows"].sum(),
                }),
                include_groups=False,
            )
            rows.append({
                "scope": scope,
                "summary": "ANI_band",
                "ani_band": band_name,
                "n_species": len(bsub),
                "n_genera": int(bg.shape[0]),
                "n_families": int(bf.shape[0]),
                "pairs": int(bsub["structural_rows"].sum()),
                "unweighted_median_fraction_ge2": float(bsub["fraction_ge2_reported"].median()),
                "pair_weighted_fraction_ge2": float(bsub["sv_ge2_pairs"].sum() / bsub["structural_rows"].sum()),
                "median_genus_weighted_fraction": float(bg["weighted_fraction"].median()),
                "median_family_weighted_fraction": float(bf["weighted_fraction"].median()),
            })
    pd.DataFrame(rows).to_csv(OUT / "gtdb_cluster_scope_sensitivity.tsv", sep="\t", index=False, float_format="%.6g")


def oapgc_alias() -> None:
    src = ROOT / "results/census/oapgc/site_effects/species_effects_all_pairs.tsv"
    df = pd.read_csv(src, sep="\t")
    df = df.rename(columns={
        "SNP_p": "ANI_distance_p",
        "SNP_cliffs_delta": "ANI_distance_cliffs_delta",
        "SNP_auc": "ANI_distance_auc",
        "SNP_q": "ANI_distance_q",
    })
    df.to_csv(ROOT / "results/census/oapgc/site_effects/species_effects_all_pairs_ani_distance.tsv", sep="\t", index=False)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    af_benchmark()
    gtdb_scope()
    oapgc_alias()


if __name__ == "__main__":
    main()
