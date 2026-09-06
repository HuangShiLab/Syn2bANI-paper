#!/usr/bin/env python3
"""Generate Supplementary Figure S12: circular-origin filtering in H. pylori."""
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SENS = ROOT / "case_studies" / "h_pylori_cagpai" / "results" / "cagpai_artifact_threshold_sensitivity.tsv"
OUT = ROOT / "paper" / "figures" / "supplementary" / "fig_s12_circular_origin_filtering.png"

before = {"empty": 85, "partial": 11, "complete_collinear": 12, "complete_rearranged": 420}
after = {"empty": 85, "partial": 11, "complete_collinear": 145, "complete_rearranged": 287}
states = list(before.keys())

sens = pd.read_csv(SENS, sep="\t")

fig, axes = plt.subplots(2, 2, figsize=(10, 9))

# Top-left: before filtering
ax = axes[0, 0]
ax.bar(range(len(states)), [before[s] for s in states], color=["#d62728", "#ff7f0e", "#2ca02c", "#1f77b4"])
ax.set_xticks(range(len(states)))
ax.set_xticklabels(states, rotation=30, ha="right")
ax.set_ylabel("Number of genomes")
ax.set_title("Before filtering")
for i, s in enumerate(states):
    ax.text(i, before[s] + 5, str(before[s]), ha="center")

# Top-right: after filtering
ax = axes[0, 1]
ax.bar(range(len(states)), [after[s] for s in states], color=["#d62728", "#ff7f0e", "#2ca02c", "#1f77b4"])
ax.set_xticks(range(len(states)))
ax.set_xticklabels(states, rotation=30, ha="right")
ax.set_ylabel("Number of genomes")
ax.set_title("After filtering")
for i, s in enumerate(states):
    ax.text(i, after[s] + 5, str(after[s]), ha="center")

# Bottom-left: rearrangement count vs threshold
ax = axes[1, 0]
ax.plot(sens["threshold"], sens["n_complete_rearranged"], "o-", color="#1f77b4", label="complete_rearranged")
ax.plot(sens["threshold"], sens["n_excluded_calls"], "s--", color="#ff7f0e", label="excluded calls")
ax.set_xlabel("Artifact span threshold")
ax.set_ylabel("Count")
ax.set_title("Rearrangement count vs. artifact threshold")
ax.legend()

# Bottom-right: p-value sensitivity
ax = axes[1, 1]
import math
pvals = sens["rearr_gc_crude_p"].replace(0, 1e-300)
neglogp = [-math.log10(p) for p in pvals]
ax.plot(sens["threshold"], neglogp, "o-", color="#1f77b4", label="rearrangement GC vs NAG")
ax.axhline(-math.log10(0.05), color="red", linestyle="--", linewidth=0.8)
ax.set_xlabel("Artifact span threshold")
ax.set_ylabel("-log10(crude p)")
ax.set_title("Association p-value sensitivity")
ax.legend()

plt.suptitle("Supplementary Figure S12 | Circular-origin artifact filtering in H. pylori cagPAI (n = 528)")
plt.tight_layout()
plt.savefig(OUT, dpi=300, bbox_inches="tight")
plt.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight")
print(f"Wrote {OUT}")
print(f"Wrote {OUT.with_suffix('.pdf')}")
