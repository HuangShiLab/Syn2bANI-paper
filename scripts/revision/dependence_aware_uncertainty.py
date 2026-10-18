#!/usr/bin/env python3
"""Dependence-aware uncertainty summaries for the Nature Methods revision.

The GTDB benchmark uses genome-level resampling: whole genomes are drawn with
replacement and every pair between two drawn genomes inherits the product of
their sampling multiplicities.  This keeps all pairs from the same genome
together.  HROM and OAPGC summaries use species-level resampling as a
transparent, computationally bounded sensitivity analysis.  It does **not**
remove shared ancestry or pair-level correlation and is reported as such.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/revision"
SEED = 20261018
B = 2000


def weighted_abs_mean(values: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sum(np.abs(values) * weights) / np.sum(weights))


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sum(values * weights) / np.sum(weights))


def genome_bootstrap(
    df: pd.DataFrame,
    estimates: dict[str, str],
    query_col: str,
    ref_col: str,
    truth_col: str,
    comparisons: list[tuple[str, str]],
) -> dict:
    genomes = pd.unique(pd.concat([df[query_col], df[ref_col]], ignore_index=True))
    pos = {g: i for i, g in enumerate(genomes)}
    gi = df[query_col].map(pos).to_numpy(np.int64)
    gj = df[ref_col].map(pos).to_numpy(np.int64)
    truth = df[truth_col].to_numpy(float)
    n = len(genomes)
    rng = np.random.default_rng(SEED)
    draws = np.empty((B, len(comparisons)), float)
    point = {}
    rows = []

    # Point estimates are restricted to the common calibrated comparisons.
    for a, b in comparisons:
        ok = np.isfinite(df[estimates[a]].to_numpy(float)) & np.isfinite(
            df[estimates[b]].to_numpy(float)
        )
        ea = df.loc[ok, estimates[a]].to_numpy(float) - truth[ok]
        eb = df.loc[ok, estimates[b]].to_numpy(float) - truth[ok]
        point[f"mae_{a}"] = float(np.mean(np.abs(ea)))
        point[f"mae_{b}"] = float(np.mean(np.abs(eb)))
        point[f"mae_{a}_minus_{b}"] = float(np.mean(np.abs(ea)) - np.mean(np.abs(eb)))

    for b in range(B):
        counts = np.bincount(rng.integers(0, n, size=n), minlength=n)
        weights = counts[gi] * counts[gj].astype(float)
        keep = weights > 0
        for j, (a, c) in enumerate(comparisons):
            ok = keep & np.isfinite(df[estimates[a]].to_numpy(float)) & np.isfinite(
                df[estimates[c]].to_numpy(float)
            )
            wa = weights[ok]
            ea = df[estimates[a]].to_numpy(float)[ok] - truth[ok]
            ec = df[estimates[c]].to_numpy(float)[ok] - truth[ok]
            ma = weighted_abs_mean(ea, wa)
            mc = weighted_abs_mean(ec, wa)
            draws[b, j] = ma - mc
    for j, (a, c) in enumerate(comparisons):
        lo, hi = np.quantile(draws[:, j], [0.025, 0.975])
        key = f"mae_{a}_minus_{c}"
        rows.append({
            "analysis": "gtdb_43k_genome_bootstrap",
            "metric": key,
            "point_estimate": point[key],
            "ci95_lower": float(lo),
            "ci95_upper": float(hi),
            "units": "ANI points; positive favours second method",
        })
    return {"point": point, "rows": rows}


def species_bootstrap(
    df: pd.DataFrame,
    value_cols: list[str],
    labels: dict[str, str],
    threshold: float = 0.5,
) -> list[dict]:
    vals = {c: df[c].to_numpy(float) for c in value_cols}
    n = len(df)
    rng = np.random.default_rng(SEED + 1)
    med = {c: np.empty(B) for c in value_cols}
    prop = {c: np.empty(B) for c in value_cols}
    for b in range(B):
        idx = rng.integers(0, n, size=n)
        counts = np.bincount(idx, minlength=n).astype(float)
        for c in value_cols:
            x = vals[c][idx]
            med[c][b] = np.median(x)
            prop[c][b] = np.sum(counts * (x > threshold)) / np.sum(counts)
    rows = []
    for c in value_cols:
        point_med = float(np.median(vals[c]))
        point_prop = float(np.mean(vals[c] > threshold))
        rows.extend([
            {
                "analysis": labels["analysis"],
                "metric": f"median_{c}",
                "point_estimate": point_med,
                "ci95_lower": float(np.quantile(med[c], 0.025)),
                "ci95_upper": float(np.quantile(med[c], 0.975)),
                "units": labels["units"],
            },
            {
                "analysis": labels["analysis"],
                "metric": f"fraction_{c}_above_{threshold:g}",
                "point_estimate": point_prop,
                "ci95_lower": float(np.quantile(prop[c], 0.025)),
                "ci95_upper": float(np.quantile(prop[c], 0.975)),
                "units": "species fraction",
            },
        ])
    return rows


def hrom_category_bootstrap(df: pd.DataFrame) -> list[dict]:
    mode = df["mode"].to_numpy()
    n = len(df)
    rng = np.random.default_rng(SEED + 2)
    cats = sorted(set(mode))
    samples = np.empty((B, len(cats)))
    for b in range(B):
        x = mode[rng.integers(0, n, size=n)]
        samples[b] = [(x == c).mean() for c in cats]
    rows = []
    for j, cat in enumerate(cats):
        rows.append({
            "analysis": "hrom_all_pairs_species_bootstrap",
            "metric": f"fraction_{cat}",
            "point_estimate": float(np.mean(mode == cat)),
            "ci95_lower": float(np.quantile(samples[:, j], 0.025)),
            "ci95_upper": float(np.quantile(samples[:, j], 0.975)),
            "units": "tested-species fraction",
        })
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    truth = pd.read_csv(ROOT / "results/gtdb50k/truth_50k.tsv", sep="\t")
    pred = pd.read_csv(ROOT / "results/gtdb50k/s2b_50k.tsv", sep="\t")
    pairmeta = pd.read_csv(ROOT / "results/gtdb50k/pairs_50k.tsv", sep="\t")
    pairmeta["pairid"] = pairmeta["q_acc"] + "__" + pairmeta["r_acc"]
    gtdb = truth.merge(pairmeta[["pairid", "q_acc", "r_acc", "skani_ani"]], on="pairid")
    gtdb = gtdb.merge(pred[["pairid", "ani_gated", "ani_cal"]], on="pairid")
    gtdb["hybrid"] = np.where(
        gtdb["ani_gated"] >= 98.0, gtdb["ani_gated"], gtdb["ani_cal"]
    )
    estimates = {"syn2bani_v6": "ani_cal", "syn2bani_hybrid": "hybrid", "skani": "skani_ani"}
    gt = genome_bootstrap(
        gtdb,
        estimates,
        "q_acc",
        "r_acc",
        "anim_ani",
        [("syn2bani_hybrid", "skani"), ("syn2bani_v6", "skani")],
    )

    oapgc = pd.read_csv(
        ROOT / "results/census/oapgc/site_effects/species_effects_all_pairs.tsv", sep="\t"
    )
    oa = species_bootstrap(
        oapgc,
        ["SNP_auc", "BP_auc"],
        {
            "analysis": "oapgc_species_bootstrap",
            "units": "Mann-Whitney AUC",
        },
    )
    # The source column retains the historical name SNP_auc, but `100 - skani
    # ANI` is an ANI-derived sequence distance, not a counted SNP distance.
    for r in oa:
        r["metric"] = r["metric"].replace("SNP_auc", "sequence_distance_auc")

    hrom = pd.read_csv(
        ROOT / "results/census/hrom/hypermodes/hypermode_species_enrichment.tsv", sep="\t"
    )
    hr = hrom_category_bootstrap(hrom)

    rows = gt["rows"] + oa + hr
    out_tsv = OUT / "dependence_aware_uncertainty.tsv"
    pd.DataFrame(rows).to_csv(out_tsv, sep="\t", index=False, float_format="%.6g")

    report = [
        "# Dependence-aware uncertainty (revision audit)",
        "",
        "## Scope and interpretation",
        "",
        "GTDB benchmark intervals use 2,000 genome-level bootstrap replicates: whole genomes are resampled with replacement and pairs inherit the product of endpoint multiplicities, so all observations containing a genome move together. HROM and OAPGC intervals use 2,000 species-level bootstrap replicates. Species-level resampling quantifies composition uncertainty but does not remove pair-level dependence, shared ancestry, assembly-quality correlation or database sampling structure; these remain explicit limitations.",
        "",
        "## GTDB-R207 43,334-pair benchmark",
        "",
        "MAE differences are reported on the common complete-case subset (39,903 pairs). Positive values mean the second method has lower MAE. The point estimates are reproducible from `results/gtdb50k/s2b_50k.tsv` with v6 calibration (`ani_cal`) and the operational hybrid rule (raw gated at ANI estimate ≥98%).",
        "",
    ]
    for r in gt["rows"]:
        report.append(
            f"- `{r['metric']}` = {r['point_estimate']:.3f} "
            f"(genome-bootstrap 95% CI {r['ci95_lower']:.3f} to {r['ci95_upper']:.3f})"
        )
    report.extend([
        "",
        "## OAPGC oral-airway site test",
        "",
        "The sequence-distance metric is `100 - skani ANI`; it is an ANI-derived distance rather than a base-level SNP count.",
        "",
    ])
    for r in oa:
        if r["metric"].startswith("median"):
            report.append(
                f"- `{r['metric']}` = {r['point_estimate']:.3f} "
                f"(species-bootstrap 95% CI {r['ci95_lower']:.3f} to {r['ci95_upper']:.3f})"
            )
    report.extend([
        "",
        "## HROM candidate divergence-mode fractions",
        "",
    ])
    for r in hr:
        report.append(
            f"- `{r['metric']}` = {r['point_estimate']:.3f} "
            f"(species-bootstrap 95% CI {r['ci95_lower']:.3f} to {r['ci95_upper']:.3f})"
        )
    report.extend([
        "",
        "Machine-readable intervals: `results/revision/dependence_aware_uncertainty.tsv`.",
        "",
    ])
    (OUT / "DEPENDENCE_AWARE_UNCERTAINTY_REPORT.md").write_text("\n".join(report))

    metadata = {
        "seed": SEED,
        "replicates": B,
        "gtdb_pairs_input": len(gtdb),
        "gtdb_common_calibrated_pairs": int(
            (gtdb["ani_cal"].notna() & gtdb["skani_ani"].notna()).sum()
        ),
        "gtdb_genomes": int(
            len(set(gtdb["q_acc"]) | set(gtdb["r_acc"]))
        ),
        "hrom_species": len(hrom),
        "oapgc_species": len(oapgc),
    }
    (OUT / "dependence_aware_uncertainty.meta.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
