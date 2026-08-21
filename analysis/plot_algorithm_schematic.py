#!/usr/bin/env python3
"""Figure 1 — Syn2bANI algorithm schematic (4 stages, synthetic toy data).

a) In-silico digestion: both genomes -> positioned Type IIB tags (4 enzymes).
b) Anchoring: shared tags (<=2 mismatches) become anchors; pigeonhole seeding.
c) Chaining: gap-penalized collinear DP, orientation blocks; an inverted
   block in the query demonstrates inversion detection.
d) Chain-restricted MLE: per-enzyme tag outcome counts inside chains
   (0/1/2 mismatches, miss) -> likelihood -> outputs.

Output: figures/report/fig1_algorithm_schematic.png/.pdf
"""
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

ROOT = Path(__file__).resolve().parent.parent
FIGD = ROOT / "figures/report"

ENZ_COLORS = {"BcgI": "#2166ac", "AlfI": "#4dac26", "AloI": "#e08214", "FalI": "#b2182b"}
enzymes = list(ENZ_COLORS)

rng = np.random.default_rng(3)
N_TAGS = 26
ref_pos = np.sort(rng.uniform(0.03, 0.97, N_TAGS))
ref_enz = rng.integers(0, 4, N_TAGS)

# query: same tags, with an inverted segment (indices 10..15) and some loss
keep = rng.random(N_TAGS) > 0.12  # ~12% tags lost
inv_lo, inv_hi = 10, 16
qry_pos = ref_pos.copy()
seg = qry_pos[inv_lo:inv_hi][::-1]
qry_pos[inv_lo:inv_hi] = seg
qry_pos += rng.normal(0, 0.004, N_TAGS)
qry_pos = np.clip(qry_pos, 0.02, 0.98)

fig = plt.figure(figsize=(13.5, 5.9))
gs = fig.add_gridspec(2, 4, height_ratios=[1.35, 0.65], hspace=0.30, wspace=0.18,
                      left=0.04, right=0.98, top=0.88, bottom=0.10)


def genome_bar(ax, y, label):
    ax.add_patch(Rectangle((0.02, y - 0.012), 0.96, 0.024, fc="#dddddd",
                           ec="k", lw=0.8, zorder=1))
    ax.text(-0.035, y, label, ha="right", va="center", fontsize=10)


def tags(ax, positions, enz_idx, y, h=0.05):
    for p, e in zip(positions, enz_idx):
        ax.plot([p, p], [y - h / 2, y + h / 2], color=ENZ_COLORS[enzymes[e]],
                lw=2.2, zorder=3)


# ---- panel a: digestion ----
axa = fig.add_subplot(gs[0, 0])
axa.set_xlim(-0.12, 1.02); axa.set_ylim(-0.15, 1.15); axa.axis("off")
axa.set_title("a  in-silico digestion", fontsize=11, loc="left")
genome_bar(axa, 0.75, "genome A")
genome_bar(axa, 0.25, "genome B")
tags(axa, ref_pos, ref_enz, 0.75)
tags(axa, qry_pos, ref_enz, 0.25)
for i, e in enumerate(enzymes):
    axa.plot([0.06 + i * 0.24, 0.10 + i * 0.24], [1.04, 1.04],
             color=ENZ_COLORS[e], lw=2.2)
    axa.text(0.105 + i * 0.24, 1.04, e, va="center", fontsize=8.5)
axa.text(0.5, -0.10, "positioned Type IIB tags (27–32 bp)", ha="center",
         fontsize=9, style="italic")

# ---- panel b: anchoring ----
axb = fig.add_subplot(gs[0, 1])
axb.set_xlim(-0.12, 1.02); axb.set_ylim(-0.15, 1.15); axb.axis("off")
axb.set_title("b  anchoring (≤ 2 mismatches)", fontsize=11, loc="left")
genome_bar(axb, 0.75, "")
genome_bar(axb, 0.25, "")
tags(axb, ref_pos, ref_enz, 0.75)
tags(axb, qry_pos, ref_enz, 0.25)
for i in range(N_TAGS):
    if keep[i]:
        axb.plot([ref_pos[i], qry_pos[i]], [0.71, 0.29], color="#888888",
                 lw=0.7, zorder=2)
    else:
        axb.plot(ref_pos[i], 0.80, marker="x", color="k", ms=5, zorder=4)
axb.text(0.5, -0.10, "shared tags → anchors\n(pigeonhole seeding, no scan)",
         ha="center", fontsize=9, style="italic")

# ---- panel c: chaining ----
axc = fig.add_subplot(gs[0, 2])
axc.set_xlim(-0.12, 1.02); axc.set_ylim(-0.15, 1.15); axc.axis("off")
axc.set_title("c  collinear chaining (2-pass)", fontsize=11, loc="left")
genome_bar(axc, 0.75, "")
genome_bar(axc, 0.25, "")
plus_col, minus_col = "#4393c3", "#d6604d"
for i in range(N_TAGS):
    if not keep[i]:
        continue
    inv = inv_lo <= i < inv_hi
    axc.plot([ref_pos[i], qry_pos[i]], [0.71, 0.29],
             color=minus_col if inv else plus_col,
             lw=1.6 if inv else 0.9, zorder=2,
             alpha=0.95 if inv else 0.75)
tags(axc, ref_pos, ref_enz, 0.75)
tags(axc, qry_pos, ref_enz, 0.25)
axc.plot([], [], color=plus_col, lw=2, label="+ chains")
axc.plot([], [], color=minus_col, lw=2, label="− (inverted) block")
axc.legend(frameon=False, loc="upper center", fontsize=8.5,
           bbox_to_anchor=(0.5, 1.04), ncols=2)
axc.text(0.5, -0.10, "gap-penalized DP;\nadaptive chain-break test",
         ha="center", fontsize=9, style="italic")

# ---- panel d: MLE ----
axd = fig.add_subplot(gs[0, 3])
counts = [820, 120, 30, 330]  # m0, m1, m2, miss (toy)
labels = ["0 mm", "1 mm", "2 mm", "miss"]
cols = ["#2166ac", "#67a9cf", "#d1e5f0", "#999999"]
bars = axd.bar(labels, counts, color=cols, ec="k", lw=0.6)
axd.set_title("d  chain-restricted MLE", fontsize=11, loc="left")
axd.set_ylabel("tags in chains")
xx = np.linspace(0, 3, 200)
ll = -((xx - 1.05) ** 2) * 350 + 1050
ax2 = axd.twinx()
ax2.plot(xx, ll, color="#b2182b", lw=1.8)
ax2.set_yticks([])
ax2.axvline(1.05, color="#b2182b", lw=0.8, ls=":")
ax2.text(1.12, 1055, "log L(â)", color="#b2182b", fontsize=9)
axd.set_ylim(0, 1150)
axd.tick_params(labelsize=9)

# ---- bottom row: outputs box spanning all columns ----
axo = fig.add_subplot(gs[1, :])
axo.axis("off")
axo.add_patch(Rectangle((0.03, 0.12), 0.60, 0.80, fc="#f4f4f4", ec="k",
                        lw=0.8, transform=axo.transAxes))
outputs = [
    "ANI ± standard error (gated gamma / uniform MLE)",
    "aligned fraction af_query / af_reference",
    "synteny_score, synteny_blocks, breakpoint_count",
    "SV calls: inversions, indels ≥ 1 kb, translocations",
    "reliability flag: ok / INCONSISTENT / BELOW_DETECTION",
    "ani_upper95 when no point estimate is responsible",
]
for i, t in enumerate(outputs):
    axo.text(0.045, 0.82 - i * 0.125, "• " + t, fontsize=10,
             transform=axo.transAxes, va="center")
arr = FancyArrowPatch((0.25, 1.10), (0.33, 0.98), transform=axo.transAxes,
                      arrowstyle="-|>", mutation_scale=18, color="k", lw=1.2)
axo.add_patch(arr)
axo.text(0.67, 0.52,
         "one sketch, one pass:\n~8 ms per pair, 58 MB peak RSS\n"
         "(vs 8–10 s per pair for dnadiff)",
         fontsize=10, transform=axo.transAxes, va="center", style="italic")

fig.suptitle("Syn2bANI: chain-restricted maximum-likelihood ANI on fixed restriction-site anchors",
             fontsize=12.5)
for ext in ["png", "pdf"]:
    fig.savefig(FIGD / f"fig1_algorithm_schematic.{ext}", dpi=300)
print("figure ->", FIGD / "fig1_algorithm_schematic.png")
