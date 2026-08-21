#!/usr/bin/env python3
"""MAG validation analysis: syn2bani / skani / FastANI vs dnadiff (ANIm) truth
on 695 MetaBAT2 bins from 35 CAMI2 samples (25 strain-madness + 10 marine).

Inputs (results/mag_validation/collect/):
  ani_fast_tools.tsv  per-pair tool output (anchor + non-anchor rows)
  truth_dnadiff.tsv   dnadiff ANIm truth for anchor pairs
  bins.tsv            per-bin assembly / CheckM2 / contamination metadata
  rep_map.tsv         GTDB representative mapping (not used for accuracy)

Outputs:
  results/mag_validation/MAG_METRICS.tsv        metric table
  figures/report/mag_validation.png/.pdf        3-panel summary figure
  stdout summary
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent.parent
COLL = ROOT / "results/mag_validation/collect"
FIGD = ROOT / "figures/report"
FIGD.mkdir(parents=True, exist_ok=True)

tools = pd.read_csv(COLL / "ani_fast_tools.tsv", sep="\t")
truth = pd.read_csv(COLL / "truth_dnadiff.tsv", sep="\t")
bins = pd.read_csv(COLL / "bins.tsv", sep="\t")

# Anchor pairs only, joined with truth and bin metadata.
df = (tools[tools.role == "anchor"]
      .merge(truth, on="bin", how="inner")
      .merge(bins, on="bin", how="left"))
df["dataset"] = df["dataset"].fillna(df["bin"].str.split("__").str[0])
for c in ["s2b_ani", "s2b_ani_gated", "skani_ani", "fastani_ani", "anim_ani",
          "s2b_ani_upper95", "n50", "checkm2_cont", "contam_pct"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

df["err_s2b"] = df["s2b_ani"] - df["anim_ani"]
df["err_s2b_gated"] = df["s2b_ani_gated"] - df["anim_ani"]
df["err_skani"] = df["skani_ani"] - df["anim_ani"]
df["err_fastani"] = df["fastani_ani"] - df["anim_ani"]
df["flagged"] = df["s2b_flag"] != "ok"

TOOLCOLS = {"syn2bani": "err_s2b", "syn2bani_raw_gated": "err_s2b_gated",
            "skani": "err_skani", "FastANI": "err_fastani"}


def metric_rows(frame, group_label, group_val):
    rows = []
    for tool, col in TOOLCOLS.items():
        e = frame[col].dropna()
        if len(e) == 0:
            continue
        rows.append({
            "group": group_label, "value": group_val, "tool": tool,
            "n": len(e),
            "MAE": e.abs().mean(),
            "median_abs_err": e.abs().median(),
            "bias": e.mean(),
            "within_0.1": (e.abs() <= 0.1).mean(),
            "within_0.5": (e.abs() <= 0.5).mean(),
            "within_1.0": (e.abs() <= 1.0).mean(),
        })
    return rows


rows = metric_rows(df, "overall", "all")
for gcol in ["dataset", "tier", "class", "af_tier"]:
    for val, sub in df.groupby(gcol):
        rows += metric_rows(sub, gcol, val)
met = pd.DataFrame(rows)
met.to_csv(ROOT / "results/mag_validation/MAG_METRICS.tsv", sep="\t",
           index=False, float_format="%.4f")

pd.set_option("display.width", 160)
print("=== overall ===")
print(met[met.group == "overall"].to_string(index=False))
print("\n=== by tier (syn2bani) ===")
print(met[(met.group == "tier") & (met.tool == "syn2bani")].to_string(index=False))
print("\n=== by class ===")
print(met[met.group == "class"].to_string(index=False))
print("\n=== by af_tier (syn2bani) ===")
print(met[(met.group == "af_tier") & (met.tool == "syn2bani")].to_string(index=False))

# --- fragmentation: |error| vs N50 / n_contigs ---
rho, p = spearmanr(df["n50"], df["err_s2b"].abs(), nan_policy="omit")
print(f"\n|s2b err| vs N50: spearman rho={rho:.3f} (p={p:.2e})")
rho2, p2 = spearmanr(df["n_contigs"], df["err_s2b"].abs(), nan_policy="omit")
print(f"|s2b err| vs n_contigs: spearman rho={rho2:.3f} (p={p2:.2e})")

# --- flag as contamination / failure detector ---
print("\n=== flag vs contamination class ===")
print(pd.crosstab(df["s2b_flag"], df["class"], margins=True))
print("\n=== flag vs |err|>0.5 ===")
big = df["err_s2b"].abs() > 0.5
print(pd.crosstab(df["flagged"], big, margins=True,
                  rownames=["flagged"], colnames=["|err|>0.5"]))
if big.sum() > 0:
    print(f"recall on |err|>0.5: {df['flagged'][big].mean():.3f}  "
          f"precision: {big[df['flagged']].mean():.3f}" if df["flagged"].sum() else "")

# --- ani_upper95 coverage ---
cov = (df["anim_ani"] <= df["s2b_ani_upper95"]).mean()
print(f"\nani_upper95 coverage (truth <= upper95): {cov:.3f}")
gap = (df["s2b_ani_upper95"] - df["s2b_ani"]).median()
print(f"median(upper95 - ani): {gap:.3f}")

# --- figure: 3 panels ---
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))

ax = axes[0]
classes = df["class"].fillna("unknown")
palette = {"clean": "#2166ac", "strain-mixed": "#e08214",
           "cross-species": "#b2182b"}
for cls, sub in df.groupby(classes):
    ax.scatter(sub["anim_ani"], sub["s2b_ani"], s=9, alpha=0.45,
               color=palette.get(cls, "#777777"), label=cls, linewidths=0)
lim = [94.5, 100.3]
ax.plot(lim, lim, color="k", lw=0.8, ls="--")
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("ANIm truth (dnadiff, %)")
ax.set_ylabel("syn2bani estimate (%)")
ax.set_title("a  syn2bani vs ANIm truth (695 MAGs)")
ax.legend(frameon=False, markerscale=2, loc="lower right")

ax = axes[1]
order = ["syn2bani", "skani", "FastANI"]
data = [df[TOOLCOLS[t]].abs().dropna() for t in order]
bp = ax.boxplot(data, tick_labels=order, showfliers=False, widths=0.55,
                patch_artist=True, medianprops=dict(color="k"))
for patch, c in zip(bp["boxes"], ["#4393c3", "#fdb863", "#b2182b"]):
    patch.set_facecolor(c); patch.set_alpha(0.7)
ax.set_ylabel("|error| vs ANIm (pp)")
ax.set_title("b  absolute error by tool")
ax.set_yscale("log")

ax = axes[2]
tiers = [t for t in ["HQ", "MQ", "LQ"] if t in set(df["tier"])]
data = [df.loc[df.tier == t, "err_s2b"].abs().dropna() for t in tiers]
bp = ax.boxplot(data, tick_labels=tiers, showfliers=False, widths=0.55,
                patch_artist=True, medianprops=dict(color="k"))
for patch in bp["boxes"]:
    patch.set_facecolor("#92c5de"); patch.set_alpha(0.8)
ax.set_ylabel("|syn2bani error| (pp)")
ax.set_title("c  error by CheckM2 quality tier")
ax.set_yscale("log")

fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(FIGD / f"mag_validation.{ext}", dpi=300)
print(f"\nfigure -> {FIGD/'mag_validation.png'}")
