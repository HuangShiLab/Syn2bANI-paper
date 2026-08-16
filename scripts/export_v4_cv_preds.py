#!/usr/bin/env python3
"""Export per-pair band-holdout CV predictions for the deployed v4 model
(seed 0, v4a feature set) — the v4 counterpart of
results/panel_by_band/ridge_cv_preds_4e.tsv, used by figure F7c.

Output columns: query reference band anim_ani ridge_pred
(query/reference are the assembly accessions where available: the old
2,053 rows carry query_asm/ref_asm from the acc2seqid map, the 467 hi95
rows are accession-keyed already.)
"""
from pathlib import Path

import pandas as pd

import calibration_v2 as c2
import calibration_v3 as c3
import calibration_v4 as c4

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results/panel_by_band/ridge_cv_preds_v4.tsv"


def main():
    old = c3.load_gated()
    new = c4.load_hi95()
    df = pd.concat([old, new], ignore_index=True)
    d = c4.shuffled(df, c4.PRIMARY_SEED)
    preds = c2.band_holdout_cv(d, c3.V3A, "ridge")

    q = d["query_asm"] if "query_asm" in d.columns else d["query"]
    r = d["ref_asm"] if "ref_asm" in d.columns else d["reference"]
    q = q.fillna(d["query"])
    r = r.fillna(d["reference"])
    out = pd.DataFrame({"query": q, "reference": r, "band": d["band"],
                        "anim_ani": d["anim_ani"], "ridge_pred": preds})
    out = out.round(4)
    out.to_csv(OUT, sep="\t", index=False)
    mae = (out["ridge_pred"] - out["anim_ani"]).abs().mean()
    print(f"wrote {OUT} (n={len(out)}, CV MAE={mae:.4f} — "
          f"should match v4a seed-0 overall 0.8523)")


if __name__ == "__main__":
    main()
