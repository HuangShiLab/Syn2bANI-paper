#!/usr/bin/env python3
"""Generate Supplementary Figure S12 for the cagPAI pilot panel."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "cagpai_pilot"
FIG = ROOT / "figures" / "report" / "fig_s12_cagpai_pilot.png"
FIG.parent.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(DATA / "dist_all.tsv", sep="\t")

# Focus on engineered-structure pairs vs wild-type
order = ["wt", "del", "inv", "transloc", "mut1", "mut1_del", "mut1_inv", "mut1_transloc"]
label_map = {
    "wt": "WT",
    "del": "ΔcagPAI",
    "inv": "cagPAI inv",
    "transloc": "cagPAI trans",
    "mut1": "mut1",
    "mut1_del": "mut1 ΔcagPAI",
    "mut1_inv": "mut1 cagPAI inv",
    "mut1_transloc": "mut1 cagPAI trans",
}

# Pairs: wt vs each other genome (query=wt)
sub = df[(df["query"] == "wt") & (df["reference"].isin(order))].copy()
sub = sub.set_index("reference").reindex(order).reset_index()
sub["label"] = sub["reference"].map(label_map)

fig, axes = plt.subplots(1, 3, figsize=(12, 4))

# Panel A: ANI (should be ~100 for same background, ~99 for mut1 background)
axes[0].bar(range(len(sub)), sub["ani_gated"], color="steelblue")
axes[0].set_xticks(range(len(sub)))
axes[0].set_xticklabels(sub["label"], rotation=45, ha="right")
axes[0].set_ylabel("ANI (%)")
axes[0].set_ylim(97.5, 100.1)
axes[0].set_title("ANI is blind to 36 kb rearrangements")

# Panel B: breakpoint count
axes[1].bar(range(len(sub)), sub["breakpoint_count"], color="coral")
axes[1].set_xticks(range(len(sub)))
axes[1].set_xticklabels(sub["label"], rotation=45, ha="right")
axes[1].set_ylabel("Breakpoint count")
axes[1].set_title("Rearrangements increase breakpoint count")

# Panel C: synteny_score (coverage) - drops modestly with deletion
axes[2].bar(range(len(sub)), sub["synteny_score"], color="seagreen")
axes[2].set_xticks(range(len(sub)))
axes[2].set_xticklabels(sub["label"], rotation=45, ha="right")
axes[2].set_ylabel("Synteny score (chain coverage)")
axes[2].set_ylim(0.99, 1.001)
axes[2].set_title("Coverage metric is insensitive to rearrangement count")

plt.tight_layout()
plt.savefig(FIG, dpi=300)
plt.savefig(FIG.with_suffix(".pdf"), dpi=300)
print(f"Saved {FIG}")
