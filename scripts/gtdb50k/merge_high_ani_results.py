#!/usr/bin/env python3
"""Merge high-ANI tool outputs into a single TSV for analysis and v6 training.

Inputs (results/gtdb50k/):
  high_ani_pairs_final.tsv   pair metadata + split + anim_ani
  s2b_high_ani_final/*.tsv   syn2bani --verbose --calibrate outputs
  skani_high_ani_final/*.txt skani dist outputs
  fastani_high_ani_final/*.txt FastANI outputs

Output:
  high_ani_results.tsv       pair metadata + syn2bani features/estimates +
                             skani + FastANI
"""
import os
import glob
import argparse
import pandas as pd
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# In the repo the script is at scripts/gtdb50k/; on HPC it may be copied to
# results/gtdb50k/scripts/. Infer the project root accordingly.
if os.path.basename(os.path.dirname(HERE)) == "gtdb50k" and os.path.basename(HERE) == "scripts":
    # HPC copy location: results/gtdb50k/scripts/<script>
    ROOT = os.environ.get("SYN2BANI_ROOT", os.path.join(HERE, "..", "..", ".."))
else:
    ROOT = os.environ.get("SYN2BANI_ROOT", os.path.join(HERE, "..", ".."))
RES = os.path.join(ROOT, "results", "gtdb50k")

S2B_COLS = [
    "query", "reference", "ani", "ani_uniform", "af_query", "af_reference",
    "std_err", "ani_cal", "synteny_blocks", "anchor_adjacency", "breakpoint_count",
    "het_shape", "retention", "ani_from_loss", "ani_from_hist",
    "enzyme_spread", "enzyme_chi2", "per_enzyme", "n_anchors", "n_chains",
    "n_tags", "max_block_anchors", "mean_block_anchors", "flag", "ani_gated",
    "gate", "ani_upper95",
]


def parse_s2b(pairid, path):
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return None
    if df.empty:
        return None
    row = df.iloc[0].to_dict()
    rec = {"pairid": pairid}
    for c in S2B_COLS:
        rec[c] = row.get(c, np.nan)
    return rec


def parse_skani(pairid, path):
    try:
        with open(path) as fh:
            header_seen = False
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                # skani 0.3.2 output: header line starts with Ref_file
                if not header_seen:
                    if parts[0].lower() in ("ref_file", "ref_name"):
                        header_seen = True
                        continue
                    header_seen = True
                if len(parts) >= 3:
                    return {"pairid": pairid,
                            "skani_ani": float(parts[2]),
                            "skani_align_frac": float(parts[3]) if len(parts) > 3 else np.nan}
    except Exception:
        pass
    return {"pairid": pairid, "skani_ani": np.nan, "skani_align_frac": np.nan}


def parse_fastani(pairid, path):
    try:
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) >= 3:
                    return {"pairid": pairid,
                            "fastani_ani": float(parts[2]),
                            "fastani_mapped": int(parts[3]) if len(parts) > 3 else np.nan,
                            "fastani_total": int(parts[4]) if len(parts) > 4 else np.nan}
    except Exception:
        pass
    return {"pairid": pairid, "fastani_ani": np.nan,
            "fastani_mapped": np.nan, "fastani_total": np.nan}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--res", default=RES)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()
    out = args.out or os.path.join(args.res, "high_ani_results.tsv")

    pairs = pd.read_csv(os.path.join(args.res, "high_ani_pairs_final.tsv"), sep="\t")

    # syn2bani
    s2b_rows = []
    for pid in pairs["pairid"]:
        path = os.path.join(args.res, "s2b_high_ani_final", f"{pid}.tsv")
        rec = parse_s2b(pid, path)
        if rec:
            s2b_rows.append(rec)
    s2b = pd.DataFrame(s2b_rows)

    # skani
    sk_rows = []
    for pid in pairs["pairid"]:
        path = os.path.join(args.res, "skani_high_ani_final", f"{pid}.txt")
        sk_rows.append(parse_skani(pid, path))
    sk = pd.DataFrame(sk_rows)

    # fastani
    fa_rows = []
    for pid in pairs["pairid"]:
        path = os.path.join(args.res, "fastani_high_ani_final", f"{pid}.txt")
        fa_rows.append(parse_fastani(pid, path))
    fa = pd.DataFrame(fa_rows)

    df = pairs.merge(s2b, on="pairid", how="left")
    df = df.merge(sk, on="pairid", how="left")
    df = df.merge(fa, on="pairid", how="left")

    df.to_csv(out, sep="\t", index=False, float_format="%.6g")
    print(f"merged {len(df)} pairs; s2b found {len(s2b)}; skani {len(sk)}; fastani {len(fa)}")
    print(f"missing s2b: {df['ani_gated'].isna().sum()}; skani: {df['skani_ani'].isna().sum()}; fastani: {df['fastani_ani'].isna().sum()}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
