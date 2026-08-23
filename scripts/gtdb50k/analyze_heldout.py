#!/usr/bin/env python3
"""Held-out evaluation of syn2bani (raw gated + v5-calibrated) and skani
against dnadiff/ANIm truth on 43,334 GTDB-R207 pairs disjoint from the
2,541-pair v5 training set.

Inputs (results/gtdb50k/):
  truth_50k.tsv  pairid, anim_ani, anim_af_ref, anim_af_qry, ...
  s2b_50k.tsv    pairid + syn2bani verbose columns (ani_gated, ani_cal, flag, ...)
  pairs_50k.tsv  q_acc, r_acc, skani_ani, band, phylum

Outputs:
  gtdb50k_metrics.tsv        overall + per-band + per-phylum metrics
  GTDB50K_HELDOUT_REPORT.md  human-readable report
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "..", "results", "gtdb50k")

BAND_ORDER = ["80-85", "85-90", "90-95", "95-100"]
METHODS = [("syn2bani_raw", "ani_gated"), ("syn2bani_cal", "ani_cal"), ("skani", "skani_ani"), ("fastani", "fastani_ani")]


def metrics(err):
    err = err[np.isfinite(err)]
    n = len(err)
    if n == 0:
        return dict(n=0, mae=np.nan, med=np.nan, bias=np.nan, sd=np.nan, r=np.nan)
    return dict(n=n, mae=np.mean(np.abs(err)), med=np.median(np.abs(err)),
                bias=np.mean(err), sd=np.std(err), r=np.nan)


def main():
    truth = pd.read_csv(os.path.join(RES, "truth_50k.tsv"), sep="\t")
    s2b = pd.read_csv(os.path.join(RES, "s2b_50k.tsv"), sep="\t")
    pairs = pd.read_csv(os.path.join(RES, "pairs_50k.tsv"), sep="\t")
    pairs["pairid"] = pairs["q_acc"] + "__" + pairs["r_acc"]

    df = truth.merge(pairs[["pairid", "skani_ani", "band", "phylum"]], on="pairid", how="inner")
    df = df.merge(s2b[["pairid", "ani_gated", "ani_cal", "ani", "gate", "flag",
                       "ani_upper95", "af_query", "synteny_score", "std_err"]],
                  on="pairid", how="inner")
    assert len(df) == len(truth), f"join dropped rows: {len(df)} vs {len(truth)}"

    fast = pd.read_csv(os.path.join(RES, "fastani_50k.tsv"), sep="\t")
    fast["fastani_ani"] = pd.to_numeric(fast["fastani_ani"], errors="coerce")
    df = df.merge(fast[["pairid", "fastani_ani"]], on="pairid", how="inner")

    for name, col in METHODS:
        df["err_" + name] = df[col] - df["anim_ani"]

    rows = []
    def block(sub, label, group):
        for name, _ in METHODS:
            err = sub["err_" + name].to_numpy(dtype=float)
            mask = np.isfinite(err)
            err = err[mask]
            m = metrics(err)
            if m["n"] > 2:
                t = sub.loc[mask, "anim_ani"].to_numpy(dtype=float)
                m["r"] = float(np.corrcoef(sub.loc[mask, ["err_" + name]].iloc[:, 0] + t, t)[0, 1])
            rows.append(dict(group=group, label=label, method=name, **m))

    block(df, "overall", "overall")
    for band in BAND_ORDER:
        block(df[df["band"] == band], band, "band")
    # common subset: only pairs where ani_cal is finite (excludes BELOW_DETECTION),
    # so skani / raw / cal are compared on exactly the same pairs
    dfc = df[np.isfinite(df["ani_cal"])]
    block(dfc, "overall", "common")
    for band in BAND_ORDER:
        block(dfc[dfc["band"] == band], band, "common_band")
    top_phyla = df["phylum"].value_counts()
    top_phyla = top_phyla[top_phyla >= 200].index
    for ph in sorted(top_phyla):
        block(df[df["phylum"] == ph], ph, "phylum")

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "gtdb50k_metrics.tsv"), sep="\t", index=False,
               float_format="%.4f")

    # diagnostics
    cal_over100 = int((df["ani_cal"] > 100).sum())
    cal_over100_by_band = df[df["ani_cal"] > 100].groupby("band").size().to_dict()
    flag_counts = df["flag"].value_counts().to_dict()
    gate_counts = df["gate"].value_counts().to_dict()
    # coverage: fraction with |anim - raw| <= (ani_upper95 - raw) i.e. truth <= upper95
    cover = float((df["anim_ani"] <= df["ani_upper95"]).mean())
    # af tiers
    af_tier = truth["af_tier"].value_counts().to_dict()

    def g(group, label, method):
        r = out[(out.group == group) & (out.label == label) & (out.method == method)]
        return r.iloc[0] if len(r) else None

    L = []
    L.append("# GTDB-R207 50k held-out benchmark report\n")
    L.append("43,334 same-genus pairs sampled from GTDB R207 representative genomes, "
             "disjoint from the 2,541-pair v5 calibration training set (both directions excluded). "
             "Truth: nucmer/dnadiff ANIm. Methods: syn2bani raw gated estimate, "
             "syn2bani v5-calibrated estimate (true extrapolation), skani, FastANI.\n")
    L.append("## Overall\n")
    L.append("| method | n | MAE | median | bias | SD | r |")
    L.append("|---|---|---|---|---|---|---|")
    for name, _ in METHODS:
        m = g("overall", "overall", name)
        L.append(f"| {name} | {m.n:.0f} | {m.mae:.4f} | {m.med:.4f} | {m.bias:+.4f} | {m.sd:.4f} | {m.r:.4f} |")
    L.append("\n## By ANI band (skani pre-screen band)\n")
    L.append("| band | method | n | MAE | bias | r |")
    L.append("|---|---|---|---|---|---|")
    for band in BAND_ORDER:
        for name, _ in METHODS:
            m = g("band", band, name)
            if m is None:
                continue
            L.append(f"| {band} | {name} | {m.n:.0f} | {m.mae:.4f} | {m.bias:+.4f} | {m.r:.4f} |")
    L.append("\n## Common subset (pairs with a calibrated score; BELOW_DETECTION excluded)\n")
    L.append("| band | method | n | MAE | bias | r |")
    L.append("|---|---|---|---|---|---|")
    for name, _ in METHODS:
        m = g("common", "overall", name)
        L.append(f"| all | {name} | {m.n:.0f} | {m.mae:.4f} | {m.bias:+.4f} | {m.r:.4f} |")
    for band in BAND_ORDER:
        for name, _ in METHODS:
            m = g("common_band", band, name)
            if m is None:
                continue
            L.append(f"| {band} | {name} | {m.n:.0f} | {m.mae:.4f} | {m.bias:+.4f} | {m.r:.4f} |")
    L.append("\n## By phylum (n>=200)\n")
    L.append("| phylum | method | n | MAE | bias |")
    L.append("|---|---|---|---|---|")
    for ph in sorted(top_phyla):
        for name, _ in METHODS:
            m = g("phylum", ph, name)
            L.append(f"| {ph} | {name} | {m.n:.0f} | {m.mae:.4f} | {m.bias:+.4f} |")
    L.append("\n## Diagnostics\n")
    L.append(f"- ani_cal > 100: {cal_over100} pairs ({100*cal_over100/len(df):.3f}%), by band: {cal_over100_by_band}")
    L.append(f"- flag distribution: {flag_counts}")
    L.append(f"- gate distribution: {gate_counts}")
    L.append(f"- upper95 coverage (anim_ani <= ani_upper95): {100*cover:.2f}%")
    L.append(f"- FastANI non-call rate by band (NA rows):")
    for band in BAND_ORDER:
        sub = df[df["band"] == band]
        rate = sub["fastani_ani"].isna().mean() if len(sub) else np.nan
        L.append(f"  - {band}: {100*rate:.1f}% ({sub['fastani_ani'].isna().sum()}/{len(sub)})")
    L.append(f"- AF tiers of truth set: {af_tier}")

    with open(os.path.join(RES, "GTDB50K_HELDOUT_REPORT.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
