#!/usr/bin/env python3
"""Synteny benchmark figure: inversion ladder (exact truth) + MAG specificity.

Panels:
  a) syn2bani breakpoint_count vs true breakpoints (2 per inversion) —
     exact recall at both ANI backgrounds.
  b) anchor_adjacency vs number of inversions (monotone decrease).
  c) ANI estimate invariance to rearrangement count (truth lines at 95/98).

Input: results/synteny_bench/synteny_ladder_results.tsv
Output: figures/report/fig_synteny_ladder.png/.pdf
"""
from pathlib import Path

import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent.parent
df = pd.read_csv(ROOT / "results/synteny_bench/synteny_ladder_results.tsv", sep="\t")
FIGD = ROOT / "figures/report"

colors = {0.95: "#2166ac", 0.98: "#b2182b"}

fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))

ax = axes[0]
for ani, sub in df.groupby("ani"):
    ax.scatter(sub.true_breakpoints, sub.breakpoint_count, s=45,
               color=colors[ani], label=f"ANI {ani:.2f}", zorder=3)
lim = [0, 70]
ax.plot(lim, lim, "k--", lw=0.8)
ax.set_xlabel("true breakpoints (2 per inversion)")
ax.set_ylabel("syn2bani breakpoint_count")
ax.set_title("a  breakpoint recall: exact")
ax.legend(frameon=False, loc="upper left")
ax.set_xlim(lim); ax.set_ylim(lim)

ax = axes[1]
for ani, sub in df.groupby("ani"):
    sub = sub.sort_values("n_inv")
    ax.plot(sub.n_inv, sub.anchor_adjacency, "o-", color=colors[ani],
            label=f"ANI {ani:.2f}")
ax.set_xscale("symlog", linthresh=1)
ax.set_xlabel("number of inversions")
ax.set_ylabel("anchor_adjacency")
ax.set_title("b  anchor adjacency vs rearrangement load")
ax.set_ylim(0.975, 1.001)

ax = axes[2]
for ani, sub in df.groupby("ani"):
    sub = sub.sort_values("n_inv")
    ax.axhline(ani, color=colors[ani], lw=0.8, ls=":")
    ax.plot(sub.n_inv, sub.ani_est, "o-", color=colors[ani],
            label=f"true ANI {ani:.2f}")
ax.set_xscale("symlog", linthresh=1)
ax.set_xlabel("number of inversions")
ax.set_ylabel("ANI estimate (%)")
ax.set_ylim(94.8, 98.2)
ax.set_title("c  ANI estimate is SV-invariant")
ax.legend(frameon=False, loc="lower left")

fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(FIGD / f"fig_synteny_ladder.{ext}", dpi=300)
print("figure ->", FIGD / "fig_synteny_ladder.png")
