#!/usr/bin/env python3
import pandas as pd
import numpy as np

df = pd.read_csv("/lustre1/g/aos_shihuang/data/validation_mid_ani/mid_ani_matrix_fastani.tsv", sep="\t")
print(f"Total within-genus pairs: {len(df)}")
print(f"Pairs with FastANI result: {df['fastani_ani'].notna().sum()}")
print(f"Pairs missing FastANI result: {df['fastani_ani'].isna().sum()}")
print()
print("FastANI result distribution by genus:")
for genus in sorted(df["q_genus"].unique()):
    sub = df[df["q_genus"] == genus]
    valid = sub["fastani_ani"].dropna()
    if len(valid) > 0:
        print(f"  {genus}: {len(sub)} pairs, {len(valid)} valid, range: {valid.min()*100:.2f}% - {valid.max()*100:.2f}%")
    else:
        print(f"  {genus}: {len(sub)} pairs, 0 valid")
print()
print("ANI bins (valid pairs only):")
valid = df.dropna(subset=["fastani_ani"]).copy()
valid["ani_pct"] = valid["fastani_ani"] * 100
bins = [0, 80, 85, 90, 95, 100]
valid["bin"] = pd.cut(valid["ani_pct"], bins=bins)
print(valid["bin"].value_counts().sort_index())
