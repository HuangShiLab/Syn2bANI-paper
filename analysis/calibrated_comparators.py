#!/usr/bin/env python3
"""Give skani and FastANI the same post-hoc calibration Syn2bANI gets.

The GTDB-R207 held-out benchmark compares Syn2bANI's ridge-calibrated column
against the *uncalibrated* outputs of skani and FastANI, whose error at 80–90%
ANIm is almost entirely a constant downward bias. That comparison is not
fair. This script fits the comparators the same way the Syn2bANI calibration
is validated — band-holdout, each band predicted from a model trained on the
other bands — using a linear model on each tool's own outputs, and reports
MAE on the same pairs.

Two subsets are reported: all 43,334 pairs, and the 39,903 "common" pairs on
which Syn2bANI returns a calibrated value (it declines the rest as
BELOW_DETECTION). The dropped pairs are summarised separately.

Inputs (results/gtdb50k/): s2b_50k.tsv, truth_50k.tsv, fastani_50k.tsv,
pairs_50k.tsv (skani pre-screen values).
Output: results/gtdb50k/calibrated_comparators.tsv and a Markdown summary.
"""
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results", "gtdb50k")
BANDS = [80, 85, 90, 95, 100.01]
LABELS = ["80-85", "85-90", "90-95", "95-100"]


def band_holdout_linear(d, feats, target="anim_ani"):
    """Predict each band with a least-squares linear model trained on the
    other bands. Returns a Series aligned to d.index."""
    pred = pd.Series(np.nan, index=d.index)
    ok = d[feats + [target]].notna().all(axis=1)
    for b in LABELS:
        tr = ok & (d.band != b)
        te = ok & (d.band == b)
        X = np.c_[d.loc[tr, feats].values, np.ones(tr.sum())]
        w = np.linalg.lstsq(X, d.loc[tr, target].values, rcond=None)[0]
        pred[te] = np.c_[d.loc[te, feats].values, np.ones(te.sum())] @ w
    return pred


def summarise(d, col, mask, name, subset):
    rows = []
    e = d[col] - d.anim_ani
    m = mask & e.notna()
    for b in ["all"] + LABELS:
        mm = m if b == "all" else (m & (d.band == b))
        if mm.sum() == 0:
            continue
        rows.append(dict(method=name, subset=subset, band=b, n=int(mm.sum()),
                         mae=e[mm].abs().mean(), bias=e[mm].mean(),
                         r=np.corrcoef(d.loc[mm, col], d.loc[mm, "anim_ani"])[0, 1]))
    return rows


def main():
    s = pd.read_csv(os.path.join(R, "s2b_50k.tsv"), sep="\t")
    t = pd.read_csv(os.path.join(R, "truth_50k.tsv"), sep="\t")
    f = pd.read_csv(os.path.join(R, "fastani_50k.tsv"), sep="\t")
    p = pd.read_csv(os.path.join(R, "pairs_50k.tsv"), sep="\t")
    p["pairid"] = p.q_acc + "__" + p.r_acc
    d = (s.merge(t, on="pairid").merge(f, on="pairid", how="left")
         .merge(p[["pairid", "skani_ani", "skani_af_min", "band"]], on="pairid", how="left"))
    assert len(d) == 43334, len(d)
    # skani pre-screen bands, as in the held-out report
    d["band"] = d["band"].astype(str)

    d["fastani_lin1"] = band_holdout_linear(d, ["fastani_ani"])
    d["fastani_lin3"] = band_holdout_linear(d, ["fastani_ani", "fastani_mapped", "fastani_total"])
    d["skani_lin1"] = band_holdout_linear(d, ["skani_ani"])
    d["skani_lin2"] = band_holdout_linear(d, ["skani_ani", "skani_af_min"])

    common = d.ani_cal.notna()
    everything = pd.Series(True, index=d.index)
    rows = []
    for col, name in [("ani_gated", "syn2bani raw gated"), ("ani_cal", "syn2bani calibrated v5"),
                      ("skani_ani", "skani"), ("skani_lin1", "skani + linear(ani), band-holdout"),
                      ("skani_lin2", "skani + linear(ani, AF), band-holdout"),
                      ("fastani_ani", "FastANI"), ("fastani_lin1", "FastANI + linear(ani), band-holdout"),
                      ("fastani_lin3", "FastANI + linear(ani, mapped, total), band-holdout")]:
        rows += summarise(d, col, everything, name, "all 43,334")
        rows += summarise(d, col, common, name, "common 39,903")
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(R, "calibrated_comparators.tsv"), sep="\t", index=False, float_format="%.4f")

    dropped = d[~common]
    md = ["# Comparators given the same calibration protocol",
          "",
          "Band-holdout linear recalibration of skani and FastANI on their own outputs, "
          "evaluated on the GTDB-R207 held-out set against ANIm. Produced by "
          "`analysis/calibrated_comparators.py`.",
          "",
          "## MAE (ANI points) by band, common subset (39,903 pairs with a Syn2bANI calibrated value)",
          ""]
    piv = out[out.subset == "common 39,903"].pivot(index="method", columns="band", values="mae")[["80-85", "85-90", "90-95", "95-100", "all"]]
    md.append(piv.round(3).to_markdown())
    md += ["", "## MAE by band, all 43,334 pairs (Syn2bANI calibrated is undefined on 3,431 of them)", ""]
    piv2 = out[out.subset == "all 43,334"].pivot(index="method", columns="band", values="mae")[["80-85", "85-90", "90-95", "95-100", "all"]]
    md.append(piv2.round(3).to_markdown())
    md += ["", "## Pearson r by band, common subset", ""]
    piv3 = out[out.subset == "common 39,903"].pivot(index="method", columns="band", values="r")[["80-85", "85-90", "90-95", "95-100", "all"]]
    md.append(piv3.round(3).to_markdown())
    md += ["", f"## The {len(dropped):,} pairs without a calibrated Syn2bANI value", "",
           f"Bands: {dropped.band.value_counts().to_dict()}",
           f"Syn2bANI raw gated MAE on them: {(dropped.ani_gated - dropped.anim_ani).abs().mean():.3f}",
           f"skani MAE on them: {(dropped.skani_ani - dropped.anim_ani).abs().mean():.3f}; "
           f"skani + linear: {(dropped.skani_lin1 - dropped.anim_ani).abs().mean():.3f}",
           f"FastANI MAE on them: {(dropped.fastani_ani - dropped.anim_ani).abs().mean():.3f}; "
           f"FastANI + linear: {(dropped.fastani_lin1 - dropped.anim_ani).abs().mean():.3f}",
           ""]
    with open(os.path.join(R, "CALIBRATED_COMPARATORS.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
