#!/usr/bin/env python3
"""Export per-pair band-holdout CV predictions for the v5 model (seed 0,
selected v5 feature set) — the v5 counterpart of
results/panel_by_band/ridge_cv_preds_v4.tsv, used by figure F7c.

Output columns: query reference band anim_ani ridge_pred
(query/reference are the assembly accessions where available: the old
2,053 rows carry query_asm/ref_asm from the acc2seqid map, the 467 hi95
rows are accession-keyed already.)
"""
import json
from pathlib import Path

import pandas as pd

import calibration_v2 as c2
import calibration_v3 as c3
import calibration_v4 as c4
import calibration_v5 as c5

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results/panel_by_band/ridge_cv_preds_v5.tsv"
V5_JSON = ROOT / "results/panel_by_band/linear_cal_v5.json"


def main():
    # feature set actually selected by the v5 run (v5a 9 feats or v5b +gate)
    n_feats = len(json.load(open(V5_JSON))["feature_names"])
    features = c3.V3A if n_feats == 9 else c3.V3B

    old = c5.load_gated_v5()
    new = c5.load_hi95_v5()
    df = pd.concat([old, new], ignore_index=True)
    d = c4.shuffled(df, c4.PRIMARY_SEED)
    preds = c2.band_holdout_cv(d, features, "ridge")

    q = d["query_asm"] if "query_asm" in d.columns else d["query"]
    r = d["ref_asm"] if "ref_asm" in d.columns else d["reference"]
    q = q.fillna(d["query"])
    r = r.fillna(d["reference"])
    out = pd.DataFrame({"query": q, "reference": r, "band": d["band"],
                        "anim_ani": d["anim_ani"], "ridge_pred": preds})
    out = out.round(4)
    out.to_csv(OUT, sep="\t", index=False)
    mae = (out["ridge_pred"] - out["anim_ani"]).abs().mean()
    print(f"wrote {OUT} (n={len(out)}, {n_feats}-feature set, "
          f"CV MAE={mae:.4f} — should match the v5 seed-0 overall)")


if __name__ == "__main__":
    main()
