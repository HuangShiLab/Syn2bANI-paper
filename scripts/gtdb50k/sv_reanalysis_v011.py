#!/usr/bin/env python3
"""SV re-analysis with the FIXED Syn2bANI v0.1.1 breakpoint_count.

Replaces the withdrawn statistics in results/gtdb50k/SV_REANALYSIS.md (computed
with the pre-v0.1.1 counter that over-counted by 1-2 orders of magnitude).

Inputs:
  results/gtdb50k/rerun_v011/s2b_50k_v011.tsv    43,334 rows, v0.1.1 binary
  results/gtdb50k/rerun_v011/s2b_pairids_sorted.txt  pairids in the merge's
      sorted-file order (1 row per file -> exact reconstruction)
  results/gtdb50k/truth_50k.tsv                  per-pair ANIm
  results/gtdb50k/dnadiff_events_50k.tsv         dnadiff event counts
  ../Syn2bANI-paper/results/gtdb50k/sv_truth_50k_min10000.tsv  size-filtered
      dnadiff breakpoints / large indels
  data/gtdb_metadata/*_contig_count_r207.tsv     assembly contig counts

Outputs results/gtdb50k/sv_reanalysis_v011.{tsv,md}
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression

ROOT = Path(__file__).resolve().parents[2]


def pcorr(d, x, y, controls):
    X = d[controls].values
    m = np.isfinite(d[x]) & np.isfinite(d[y]) & np.all(np.isfinite(X), axis=1)
    if m.sum() < 10:
        return np.nan
    Xc = X[m]
    xr = d.loc[m, x].values - LinearRegression().fit(Xc, d.loc[m, x].values).predict(Xc)
    yr = d.loc[m, y].values - LinearRegression().fit(Xc, d.loc[m, y].values).predict(Xc)
    return pearsonr(xr, yr)[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--min10000",
                    default=str(ROOT.parent / "Syn2bANI-paper" / "results" /
                                "gtdb50k" / "sv_truth_50k_min10000.tsv"))
    args = ap.parse_args()
    res = ROOT / "results" / "gtdb50k"

    # ---- attach pairids to the merged v0.1.1 table ------------------------
    s2b = pd.read_csv(res / "rerun_v011" / "s2b_50k_v011.tsv", sep="\t")
    pairids = [l.strip() for l in
               open(res / "rerun_v011" / "s2b_pairids_sorted.txt") if l.strip()]
    assert len(s2b) == len(pairids), (len(s2b), len(pairids))
    s2b["pairid"] = pairids

    truth = pd.read_csv(res / "truth_50k.tsv", sep="\t")[["pairid", "anim_ani"]]
    ev = pd.read_csv(res / "dnadiff_events_50k.tsv", sep="\t")
    l10k = pd.read_csv(args.min10000, sep="\t")
    df = (s2b.merge(truth, on="pairid", how="left")
             .merge(ev, on="pairid", how="left")
             .merge(l10k, on="pairid", how="left"))
    df["dd_bp"] = df[["dd_breakpoints_ref", "dd_breakpoints_qry"]].max(axis=1)
    df["dd_inv"] = df[["dd_inversions_ref", "dd_inversions_qry"]].max(axis=1)

    def strip(a):
        return a[3:] if str(a).startswith(("GB_", "RS_")) else a

    cc = pd.concat([
        pd.read_csv(ROOT / "data/gtdb_metadata/ar53_contig_count_r207.tsv", sep="\t"),
        pd.read_csv(ROOT / "data/gtdb_metadata/bac120_contig_count_r207.tsv", sep="\t")
    ]).drop_duplicates("accession")
    cc["acc"] = cc.accession.map(strip)
    cmap = dict(zip(cc.acc, cc.contig_count))
    parts = df.pairid.str.split("__", expand=True)
    df["q_acc"], df["r_acc"] = parts[0], parts[1]
    df["q_contigs"] = df.q_acc.map(cmap)
    df["r_contigs"] = df.r_acc.map(cmap)
    df["n_contigs"] = df[["q_contigs", "r_contigs"]].max(axis=1)
    df.to_csv(res / "rerun_v011" / "s2b_50k_v011_merged.tsv", sep="\t", index=False)

    # sanity: any self pairs / zero checks
    n_self = int((df.q_acc == df.r_acc).sum())
    print(f"n={len(df):,}  self-pairs={n_self}  "
          f"breakpoint_count: median={df.breakpoint_count.median():.0f} "
          f"mean={df.breakpoint_count.mean():.1f} max={df.breakpoint_count.max():.0f}")

    # ---- correlations ------------------------------------------------------
    rows = []
    targets = [("dd_bp", "dnadiff breakpoints (all)"),
               ("dnadiff_breakpoints_min10000", "dnadiff breakpoints (>=10 kb)"),
               ("dnadiff_large_indels_min10000", "dnadiff large indels (>=10 kb)"),
               ("dd_inv", "dnadiff inversions (all)")]
    for col, label in targets:
        if col not in df:
            continue
        rows.append({
            "target": label,
            "raw_r": pearsonr(df.breakpoint_count, df[col])[0],
            "anim_only": pcorr(df, "breakpoint_count", col, ["anim_ani"]),
            "contigs_only": pcorr(df, "breakpoint_count", col, ["n_contigs"]),
            "both": pcorr(df, "breakpoint_count", col, ["anim_ani", "n_contigs"]),
            "n": int(df[[col]].notna().sum().iloc[0] if hasattr(df[[col]].notna().sum(), "iloc") else df[[col]].notna().sum()),
        })
    out = pd.DataFrame(rows)

    # fragmentation dependence of the two count metrics
    frag = pd.DataFrame([
        {"metric": "v0.1.1 breakpoint_count",
         "r_vs_contigs": pearsonr(df.breakpoint_count, df.n_contigs)[0],
         "spearman_vs_contigs": spearmanr(df.breakpoint_count, df.n_contigs).statistic},
        {"metric": "v0.1.1 synteny_blocks",
         "r_vs_contigs": pearsonr(df.synteny_blocks, df.n_contigs)[0],
         "spearman_vs_contigs": spearmanr(df.synteny_blocks, df.n_contigs).statistic},
    ])

    # >=95% ANIm subset (the old rho=0.674 claim)
    hi = df[df.anim_ani >= 95]
    rho_hi = spearmanr(hi.breakpoint_count, hi.dnadiff_large_indels_min10000).statistic

    # ---- cross-implementation concordance (NEW) -----------------------------
    syn2b = pd.read_csv(res / "inverted_fraction_truth_four.tsv", sep="\t")
    x = df.merge(syn2b[["pairid", "syn2b_breakpoints"]], on="pairid")
    conc = (spearmanr(x.breakpoint_count, x.syn2b_breakpoints).statistic,
            pearsonr(x.breakpoint_count, x.syn2b_breakpoints)[0], len(x))
    med_ratio = (x.breakpoint_count / x.syn2b_breakpoints.replace(0, np.nan)).median()

    # ---- write --------------------------------------------------------------
    out.to_csv(res / "sv_reanalysis_v011.tsv", sep="\t", index=False)
    lines = [
        "# SV re-analysis — Syn2bANI v0.1.1 breakpoint_count (43,334 pairs)",
        "",
        f"- breakpoint_count: median {df.breakpoint_count.median():.0f}, "
        f"mean {df.breakpoint_count.mean():.1f}, max {df.breakpoint_count.max():.0f}",
        f"- correlation with contig count: r = {frag.iloc[0].r_vs_contigs:.3f}, "
        f"rho = {frag.iloc[0].spearman_vs_contigs:.3f} (v0.1.1 counter)",
        f"- synteny_blocks vs contig count: r = {frag.iloc[1].r_vs_contigs:.3f} "
        f"(new primary-block definition)",
        "",
        "| target | raw r | \\|ANIm | \\|n_contigs | \\|both | n |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for _, r0 in out.iterrows():
        lines.append(f"| {r0.target} | {r0.raw_r:.3f} | {r0.anim_only:.3f} | "
                     f"{r0.contigs_only:.3f} | **{r0.both:.3f}** | {r0.n:,} |")
    lines += [
        "",
        f"- >=95% ANIm subset (n = {len(hi):,}): Spearman rho vs large indels = "
        f"{rho_hi:.3f} (withdrawn value was 0.674)",
        f"- cross-implementation concordance on the same {conc[2]:,} pairs: "
        f"Spearman rho = {conc[0]:.3f}, Pearson r = {conc[1]:.3f}; "
        f"median syn2bani/syn2b count ratio = {med_ratio:.2f}",
        "",
        "Replaces the withdrawn r = 0.414/0.453 table in SV_REANALYSIS.md.",
    ]
    (res / "sv_reanalysis_v011.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
