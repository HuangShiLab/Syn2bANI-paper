#!/usr/bin/env python3
"""Unified GTDB-R207 benchmark: combine 43,334 held-out (80-95%) and 2,342
high-ANI (95-100%) pairs, compute metrics, and generate figures.

Outputs:
  results/gtdb50k/unified_metrics.tsv
  results/gtdb50k/unified_benchmark_report.md
  figures/gtdb_r207_unified_scatter.png
  figures/gtdb_r207_unified_error_by_band.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("SYN2BANI_ROOT", os.path.join(HERE, "..", ".."))
RES = os.path.join(ROOT, "results", "gtdb50k")
FIG = os.path.join(ROOT, "figures")

BAND_ORDER = ["80-85", "85-90", "90-95", "95-97", "97-100"]
METHODS = {
    "syn2bani_raw": "ani_gated",
    "syn2bani_v6": "ani_cal",
    "skani": "skani_ani",
    "fastani": "fastani_ani",
}


def metrics(err):
    err = err[np.isfinite(err)]
    if len(err) == 0:
        return dict(n=0, mae=np.nan, med=np.nan, bias=np.nan, sd=np.nan, r=np.nan)
    return dict(n=len(err), mae=np.mean(np.abs(err)), med=np.median(np.abs(err)),
                bias=np.mean(err), sd=np.std(err),
                r=np.nan)


def load_heldout():
    truth = pd.read_csv(os.path.join(RES, "truth_50k.tsv"), sep="\t")
    s2b = pd.read_csv(os.path.join(RES, "s2b_50k.tsv"), sep="\t")
    pairs = pd.read_csv(os.path.join(RES, "pairs_50k.tsv"), sep="\t")
    pairs["pairid"] = pairs["q_acc"] + "__" + pairs["r_acc"]
    df = truth.merge(pairs[["pairid", "skani_ani", "band", "phylum"]], on="pairid", how="inner")
    df = df.merge(s2b[["pairid", "ani_gated", "ani_cal", "flag", "gate"]], on="pairid", how="inner")
    fast = pd.read_csv(os.path.join(RES, "fastani_50k.tsv"), sep="\t")
    fast["fastani_ani"] = pd.to_numeric(fast["fastani_ani"], errors="coerce")
    df = df.merge(fast[["pairid", "fastani_ani"]], on="pairid", how="inner")
    df["source"] = "heldout_43k"
    return df


def load_high_ani():
    df = pd.read_csv(os.path.join(RES, "high_ani_results.tsv"), sep="\t")
    df = df[df["split"] == "test"].copy()
    df["source"] = "high_ani_test"
    df["phylum"] = np.nan
    return df


def main():
    hld = load_heldout()
    hai = load_high_ani()

    # unify bands: heldout already has 80-85/85-90/90-95/95-100; high-ani has 95-97/97-100
    # split heldout 95-100 into 95-97 and 97-100 by anim_ani
    def refine_band(r):
        if r["band"] == "95-100" and r["source"] == "heldout_43k":
            if r["anim_ani"] < 97.0:
                return "95-97"
            return "97-100"
        return r["band"]
    hld["band"] = hld.apply(refine_band, axis=1)

    all_cols = set(hld.columns) & set(hai.columns)
    unified = pd.concat([hld[list(all_cols)], hai[list(all_cols)]], ignore_index=True)

    rows = []
    def block(sub, label, group):
        for name, col in METHODS.items():
            if col not in sub.columns:
                continue
            err = (sub[col] - sub["anim_ani"]).to_numpy(dtype=float)
            mask = np.isfinite(err) & sub["anim_ani"].notna().to_numpy()
            err = err[mask]
            m = metrics(err)
            if m["n"] > 2:
                t = sub.loc[mask, "anim_ani"].to_numpy(dtype=float)
                p = sub.loc[mask, col].to_numpy(dtype=float)
                m["r"] = float(np.corrcoef(p, t)[0, 1])
            rows.append(dict(group=group, label=label, method=name, **m))

    block(unified, "overall", "overall")
    for band in BAND_ORDER:
        block(unified[unified["band"] == band], band, "band")
    # common subset where ani_cal is finite
    dfc = unified[np.isfinite(unified["ani_cal"])]
    block(dfc, "overall", "common")
    for band in BAND_ORDER:
        block(dfc[dfc["band"] == band], band, "common_band")

    metrics_df = pd.DataFrame(rows)
    metrics_df.to_csv(os.path.join(RES, "unified_metrics.tsv"), sep="\t", index=False, float_format="%.4f")

    # report
    def g(group, label, method):
        r = metrics_df[(metrics_df.group == group) & (metrics_df.label == label) & (metrics_df.method == method)]
        return r.iloc[0] if len(r) else None

    L = []
    L.append("# Unified GTDB-R207 benchmark report\n")
    L.append("Combined 43,334 held-out same-genus pairs (80-95% plus 95-100 split) "
             "with 727 high-ANI test pairs (95-97 / 97-100). Truth: dnadiff/ANIm.\n")
    L.append("## Overall (all pairs with a method estimate)\n")
    L.append("| method | n | MAE | median | bias | SD | r |")
    L.append("|---|---|---|---|---|---|---|")
    for name in METHODS:
        m = g("overall", "overall", name)
        if m is None:
            continue
        L.append(f"| {name} | {m.n:.0f} | {m.mae:.4f} | {m.med:.4f} | {m.bias:+.4f} | {m.sd:.4f} | {m.r:.4f} |")
    L.append("\n## By ANI band\n")
    L.append("| band | method | n | MAE | bias | r |")
    L.append("|---|---|---|---|---|---|")
    for band in BAND_ORDER:
        for name in METHODS:
            m = g("band", band, name)
            if m is None:
                continue
            L.append(f"| {band} | {name} | {m.n:.0f} | {m.mae:.4f} | {m.bias:+.4f} | {m.r:.4f} |")
    L.append("\n## Notes\n")
    L.append("- syn2bani_v6 is the ridge calibration trained on 2,520 v5 pairs + 1,614 high-ANI train pairs.\n")
    L.append("- syn2bani_raw is the gated MLE without calibration.\n")
    L.append("- skani/FastANI are shown where they returned a value.\n")
    with open(os.path.join(RES, "UNIFIED_BENCHMARK_REPORT.md"), "w") as fh:
        fh.write("\n".join(L))

    # plots
    os.makedirs(FIG, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    colors = {"syn2bani_raw": "C0", "syn2bani_v6": "C1", "skani": "C2", "fastani": "C3"}
    # scatter
    ax = axes[0, 0]
    for name, col in METHODS.items():
        sub = unified[unified[col].notna() & unified["anim_ani"].notna()]
        ax.scatter(sub["anim_ani"], sub[col], s=3, alpha=0.3, c=colors[name], label=name)
    ax.plot([80, 100], [80, 100], "k--", lw=1)
    ax.set_xlabel("ANIm truth (%)")
    ax.set_ylabel("Estimated ANI (%)")
    ax.set_xlim(80, 100)
    ax.set_ylim(80, 100)
    ax.legend(loc="lower right", markerscale=3)
    ax.set_title("Unified GTDB-R207 benchmark")
    # error by band
    ax = axes[0, 1]
    x = np.arange(len(BAND_ORDER))
    width = 0.2
    for i, name in enumerate(METHODS):
        vals = []
        for band in BAND_ORDER:
            m = g("band", band, name)
            vals.append(m.mae if m is not None else np.nan)
        ax.bar(x + i * width, vals, width, label=name, color=colors[name])
    ax.set_xticks(x + width * 1.5)
    ax.set_xticklabels(BAND_ORDER)
    ax.set_ylabel("MAE (ANI points)")
    ax.set_xlabel("ANI band")
    ax.legend()
    ax.set_title("MAE by band")
    # signed error distribution
    ax = axes[1, 0]
    for name, col in METHODS.items():
        err = (unified[col] - unified["anim_ani"]).dropna()
        ax.hist(err, bins=50, alpha=0.4, label=name, color=colors[name])
    ax.axvline(0, color="k", linestyle="--")
    ax.set_xlabel("Error (estimate - truth, ANI points)")
    ax.set_ylabel("Count")
    ax.legend()
    ax.set_title("Signed-error distribution")
    # high-ANI zoom
    ax = axes[1, 1]
    sub = unified[unified["band"].isin(["95-97", "97-100"])]
    for name, col in METHODS.items():
        s = sub[sub[col].notna() & sub["anim_ani"].notna()]
        ax.scatter(s["anim_ani"], s[col], s=5, alpha=0.4, c=colors[name], label=name)
    ax.plot([95, 100], [95, 100], "k--", lw=1)
    ax.set_xlabel("ANIm truth (%)")
    ax.set_ylabel("Estimated ANI (%)")
    ax.set_xlim(95, 100)
    ax.set_ylim(95, 100)
    ax.legend(loc="lower right")
    ax.set_title("High-ANI regime (95-100%)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "gtdb_r207_unified_benchmark.png"), dpi=300)
    plt.close(fig)
    print("wrote", os.path.join(RES, "unified_metrics.tsv"))
    print("wrote", os.path.join(RES, "UNIFIED_BENCHMARK_REPORT.md"))
    print("wrote", os.path.join(FIG, "gtdb_r207_unified_benchmark.png"))


if __name__ == "__main__":
    main()
