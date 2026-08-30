#!/usr/bin/env python3
"""Re-analysis of the GTDB50k SV comparison, controlling for the two confounders
that the first pass did not separate: sequence divergence and assembly
fragmentation.

The first analysis (SV_COMPARISON_REPORT.md, SV_DNADIFF_FILTERED_CORRELATION.md)
reported raw correlations between the Syn2bANI structural columns and dnadiff.
Raw correlation cannot distinguish "both methods measure rearrangement" from
"both methods measure divergence" or "both methods measure contig count", and at
85-95% ANI on draft assemblies all three are in play.

Contig count is not in the merged table, but it is recoverable exactly:
breakpoint_count is defined as n_chains - n_chained_contigs and synteny_blocks is
n_chains, so

    n_chained_contigs = synteny_blocks - breakpoint_count

An independent divergence axis comes from truth_50k.tsv (ANIm), rather than from
Syn2bANI's own ani_gated, which would make the control partly circular.

Usage:  python3 scripts/sv_reanalysis.py [results/gtdb50k]
"""
import sys
import numpy as np
import pandas as pd
from scipy import stats

D = sys.argv[1] if len(sys.argv) > 1 else "results/gtdb50k"


def load():
    d = pd.read_csv(f"{D}/sv_comparison_merged.tsv", sep="\t")
    for c in d.columns[1:]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    t = pd.read_csv(f"{D}/truth_50k.tsv", sep="\t")
    d = d.merge(t[["pairid", "anim_ani"]], on="pairid", how="left")
    # Chained query contigs, exact from the definition of breakpoint_count.
    d["n_ctg"] = d["synteny_blocks"] - d["breakpoint_count"]
    return d.dropna(subset=["anim_ani", "breakpoint_count", "synteny_blocks"])


def partial(x, y, ctrls):
    """Pearson/Spearman of x and y after regressing both on the controls."""
    m = ~(x.isna() | y.isna())
    for c in ctrls:
        m &= ~c.isna()
    Z = np.column_stack([np.ones(m.sum())] + [c[m].values for c in ctrls])
    rx = x[m].values - Z @ np.linalg.lstsq(Z, x[m].values, rcond=None)[0]
    ry = y[m].values - Z @ np.linalg.lstsq(Z, y[m].values, rcond=None)[0]
    return stats.pearsonr(rx, ry)[0], stats.spearmanr(rx, ry)[0], int(m.sum())


def main():
    d = load()
    rows = []
    print(f"n = {len(d)}\n")

    print("=== 1. ANIm composition of the validation set ===")
    for lo, hi in [(0, 85), (85, 90), (90, 95), (95, 97), (97, 99), (99, 101)]:
        s = d[(d.anim_ani >= lo) & (d.anim_ani < hi)]
        print(f"  ANIm {lo:>3}-{hi:<3}: n={len(s):>6}  ({100*len(s)/len(d):>5.1f}%)")
    print(f"  >= 97% ANIm (the strain range the tool targets): n={(d.anim_ani>=97).sum()}\n")

    print("=== 2. What each column tracks ===")
    print(f"{'column':>30} {'r vs ANIm':>10} {'r vs n_ctg':>11} {'rho vs n_ctg':>13}")
    cols = ["breakpoint_count", "synteny_blocks", "dnadiff_breakpoints",
            "dnadiff_breakpoints_min10000", "dnadiff_blocks",
            "dnadiff_large_indels_min10000", "mm2_breakpoints",
            "anchor_adjacency", "af_query"]
    for c in cols:
        if c not in d:
            continue
        m = ~d[c].isna()
        ra = stats.pearsonr(d.loc[m, "anim_ani"], d.loc[m, c])[0]
        rc = stats.pearsonr(d.loc[m, "n_ctg"], d.loc[m, c])[0]
        oc = stats.spearmanr(d.loc[m, "n_ctg"], d.loc[m, c])[0]
        print(f"{c:>30} {ra:>10.3f} {rc:>11.3f} {oc:>13.3f}")
        rows.append(("tracks", c, "", ra, rc, oc, int(m.sum())))
    print(f"\n  contig count: mean {d.n_ctg.mean():.1f}, median {d.n_ctg.median():.0f}, "
          f"max {d.n_ctg.max():.0f}")
    print(f"  synteny_blocks mean {d.synteny_blocks.mean():.1f} -> "
          f"{100*d.n_ctg.mean()/d.synteny_blocks.mean():.0f}% of 'synteny blocks' are contig starts\n")

    print("=== 3. Does the agreement survive the controls? ===")
    print(f"{'pair':>54} {'raw':>7} {'|ANIm':>7} {'|n_ctg':>7} {'|both':>7}")
    pairs = [("breakpoint_count", "dnadiff_breakpoints"),
             ("breakpoint_count", "dnadiff_breakpoints_min10000"),
             ("breakpoint_count", "dnadiff_large_indels_min10000"),
             ("synteny_blocks", "dnadiff_blocks"),
             ("synteny_blocks", "dnadiff_breakpoints"),
             ("breakpoint_count", "mm2_breakpoints")]
    for a, b in pairs:
        if b not in d:
            continue
        m = ~(d[a].isna() | d[b].isna())
        raw = stats.pearsonr(d.loc[m, a], d.loc[m, b])[0]
        r1, _, _ = partial(d[a], d[b], [d.anim_ani])
        r2, _, _ = partial(d[a], d[b], [d.n_ctg])
        r3, o3, n = partial(d[a], d[b], [d.anim_ani, d.n_ctg])
        print(f"{a+' ~ '+b:>54} {raw:>7.3f} {r1:>7.3f} {r2:>7.3f} {r3:>7.3f}")
        rows.append(("partial", a, b, raw, r1, r2, r3))

    print("\n=== 4. Is the truth axis clean? ===")
    sl, ic, r, _, _ = stats.linregress(d.breakpoint_count, d.dnadiff_breakpoints)
    print(f"  dnadiff_breakpoints = {sl:.2f} x breakpoint_count + {ic:.1f}  (r={r:.3f})")
    z = d[d.breakpoint_count == 0]
    print(f"  where breakpoint_count == 0 (n={len(z)}): dnadiff_breakpoints median "
          f"{z.dnadiff_breakpoints.median():.0f}, mean {z.dnadiff_breakpoints.mean():.0f}, "
          f"at median ANIm {z.anim_ani.median():.2f}")
    m = ~d.dnadiff_breakpoints.isna()
    X = np.column_stack([np.ones(m.sum()), d.loc[m, "anim_ani"], d.loc[m, "n_ctg"],
                         d.loc[m, "dnadiff_blocks"]])
    y = d.loc[m, "dnadiff_breakpoints"].values
    pred = X @ np.linalg.lstsq(X, y, rcond=None)[0]
    r2_full = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    r2_blocks = stats.pearsonr(d.loc[m, "dnadiff_blocks"], y)[0] ** 2
    print(f"  R^2 of dnadiff_breakpoints ~ ANIm + n_ctg + dnadiff_blocks: {r2_full:.3f}")
    print(f"  R^2 from dnadiff_blocks alone:                              {r2_blocks:.3f}")

    print("\n=== 5. The >= 95% ANIm subset, where rearrangement should dominate ===")
    h = d[d.anim_ani >= 95]
    print(f"  n = {len(h)}")
    for a, b in [("breakpoint_count", "dnadiff_breakpoints"),
                 ("breakpoint_count", "dnadiff_large_indels_min10000"),
                 ("synteny_blocks", "dnadiff_blocks")]:
        mm = ~(h[a].isna() | h[b].isna())
        print(f"  {a} ~ {b}: r={stats.pearsonr(h.loc[mm,a],h.loc[mm,b])[0]:.3f} "
              f"rho={stats.spearmanr(h.loc[mm,a],h.loc[mm,b])[0]:.3f}")
    print(f"  breakpoint_count vs n_ctg here: r={stats.pearsonr(h.breakpoint_count,h.n_ctg)[0]:.3f}")
    print(f"  synteny_blocks   vs n_ctg here: r={stats.pearsonr(h.synteny_blocks,h.n_ctg)[0]:.3f}")

    out = pd.DataFrame(rows, columns=["kind", "x", "y", "v1", "v2", "v3", "v4"])
    out.to_csv(f"{D}/sv_reanalysis_metrics.tsv", sep="\t", index=False)
    print(f"\nwrote {D}/sv_reanalysis_metrics.tsv")


main()
