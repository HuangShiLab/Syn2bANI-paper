#!/usr/bin/env python3
"""Build the three joined evaluation tables used by the gating/flag prototypes.

Outputs (written next to this script):
  gtdb_anim_joined.tsv   - v8current features + ANIm truth + band (accession-keyed join)
  oral_gut_joined.tsv    - v8current features + FastANI reference (100 same-species pairs)
  midani_joined.tsv      - 15 mid-ANI pairs with ANIm truth (already keyed by accession)
"""
import pathlib

import pandas as pd

RES = pathlib.Path(__file__).resolve().parent.parent
DATA = RES.parent / "data"
OUT = pathlib.Path(__file__).resolve().parent


def gtdb():
    feat = pd.read_csv(RES / "anim_truth_2074_v8current.tsv", sep="\t")
    m = pd.read_csv(RES / "anim_2074_acc2seqid.tsv", sep="\t", header=None,
                    names=["accession", "seqid"])
    ev = pd.read_csv(RES / "panel_by_band/eval_pairs.tsv", sep="\t")

    f = feat.merge(m, left_on="query", right_on="seqid").drop(columns=["seqid"]) \
            .rename(columns={"accession": "query_asm"})
    f = f.merge(m, left_on="reference", right_on="seqid").drop(columns=["seqid"]) \
         .rename(columns={"accession": "ref_asm"})
    j = f.merge(ev[["query_asm", "ref_asm", "band", "anim_ani", "skani_ani"]],
                on=["query_asm", "ref_asm"], how="inner")
    print(f"GTDB joined: {len(j)} rows "
          f"(features {len(feat)}, finite ani {feat['ani'].notna().sum()})")
    j.to_csv(OUT / "gtdb_anim_joined.tsv", sep="\t", index=False)


def oral_gut():
    feat = pd.read_csv(RES / "oral_gut_1225_v8current.tsv", sep="\t")
    m = pd.read_csv(RES / "oral_gut_1225_acc2seqid.tsv", sep="\t", header=None,
                    names=["accession", "seqid"])
    val = pd.read_csv(DATA / "oral_gut_validation_merged_v8.tsv", sep="\t")

    f = feat.merge(m, left_on="query", right_on="seqid").drop(columns=["seqid"]) \
            .rename(columns={"accession": "query_acc"})
    f = f.merge(m, left_on="reference", right_on="seqid").drop(columns=["seqid"]) \
         .rename(columns={"accession": "ref_acc"})
    v = val[["query", "reference", "label", "fastani_ani"]].rename(
        columns={"query": "query_acc", "reference": "ref_acc"})
    j = f.merge(v, on=["query_acc", "ref_acc"], how="inner")
    # FastANI reference is on the 0-1 scale in this table
    if j["fastani_ani"].max() <= 1.5:
        j["fastani_ani"] = j["fastani_ani"] * 100.0
    print(f"oral/gut joined: {len(j)} rows; labels: {j['label'].value_counts().to_dict()}")
    j.to_csv(OUT / "oral_gut_joined.tsv", sep="\t", index=False)


def midani():
    j = pd.read_csv(RES / "validation_mid_ani_anim/anim_4e/anim_midani_evaluation.tsv", sep="\t")
    print(f"mid-ANI: {len(j)} rows")
    j.to_csv(OUT / "midani_joined.tsv", sep="\t", index=False)


if __name__ == "__main__":
    gtdb()
    oral_gut()
    midani()
