#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sens = pd.read_csv(ROOT / "results" / "cagpai_artifact_threshold_sensitivity.tsv", sep="	")

fig, axes = plt.subplots(1, 2, figsize=(9, 3.5), constrained_layout=True)

ax = axes[0]
ax.plot(sens["threshold"], sens["n_complete_rearranged"], "o-", label="complete_rearranged")
ax.plot(sens["threshold"], sens["n_excluded_calls"], "s--", label="excluded calls")
ax.set_xlabel("Artifact span threshold")
ax.set_ylabel("Count")
ax.set_title("cagPAI rearrangement count vs. artifact threshold")
ax.legend()

ax = axes[1]
import math
pvals = sens["rearr_gc_crude_p"].replace(0, 1e-300)
neglogp = [-math.log10(p) for p in pvals]
ax.plot(sens["threshold"], neglogp, "o-", label="rearrangement GC vs NAG")
ax.set_xlabel("Artifact span threshold")
ax.set_ylabel("-log10(p) crude")
ax.set_title("Association p-value sensitivity")
ax.axhline(-math.log10(0.05), color="red", linestyle="--", linewidth=0.8)
ax.legend()

fig.savefig(ROOT / "results" / "cagpai_artifact_threshold_sensitivity.png", dpi=300)
print("Wrote", ROOT / "results" / "cagpai_artifact_threshold_sensitivity.png")
