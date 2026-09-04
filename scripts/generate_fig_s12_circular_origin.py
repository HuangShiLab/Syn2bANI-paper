#!/usr/bin/env python3
"""Generate Supplementary Figure S12: circular-origin filtering in H. pylori."""
import matplotlib.pyplot as plt
import pandas as pd

before = {"empty": 85, "partial": 11, "complete_collinear": 12, "complete_rearranged": 420}
after = {"empty": 85, "partial": 11, "complete_collinear": 145, "complete_rearranged": 287}
states = list(before.keys())

fig, axes = plt.subplots(1, 2, figsize=(10, 5))

for ax, title, data in zip(axes, ["Before filtering", "After filtering"], [before, after]):
    ax.bar(range(len(states)), [data[s] for s in states], color=["#d62728", "#ff7f0e", "#2ca02c", "#1f77b4"])
    ax.set_xticks(range(len(states)))
    ax.set_xticklabels(states, rotation=30, ha="right")
    ax.set_ylabel("Number of genomes")
    ax.set_title(title)
    for i, s in enumerate(states):
        ax.text(i, data[s] + 5, str(data[s]), ha="center")

plt.suptitle("Supplementary Figure S12 | Circular-origin artifact filtering in H. pylori cagPAI (n = 528)")
plt.tight_layout()
out = "paper/figures/supplementary/fig_s12_circular_origin_filtering.png"
plt.savefig(out, dpi=300, bbox_inches="tight")
print(f"Wrote {out}")
