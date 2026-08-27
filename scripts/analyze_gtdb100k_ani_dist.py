#!/usr/bin/env python3
import pandas as pd
import numpy as np

df = pd.read_csv("/lustre1/g/aos_shihuang/Syn2bANI-paper/results/matrix_gtdb_r207_100k.tsv", sep="\t")
print(f"Total rows: {len(df)}")
print(f"Rows with FastANI: {df['fastani_ani'].notna().sum()}")
print()
print("Label distribution:")
print(df["label"].value_counts())
print()
print("FastANI ANI distribution (valid rows only):")
valid = df.dropna(subset=["fastani_ani"]).copy()
valid["ani_pct"] = valid["fastani_ani"] * 100
print(valid["ani_pct"].describe())
print()
print("ANI bins:")
bins = [0, 80, 85, 90, 95, 100]
valid["bin"] = pd.cut(valid["ani_pct"], bins=bins)
print(valid["bin"].value_counts().sort_index())
print()
print("By label, FastANI range:")
for label in sorted(valid["label"].unique()):
    sub = valid[valid["label"] == label]
    print(f"  {label}: n={len(sub)}, ANI range {sub['ani_pct'].min():.2f}% - {sub['ani_pct'].max():.2f}%, mean {sub['ani_pct'].mean():.2f}%")

print()
print("=== Error analysis by ANI bin ===")
valid["s2b_error"] = (valid["s2b_raw_ani"] - valid["fastani_ani"]) * 100
valid["s2b_v4_error"] = (valid["s2b_gbrt_ani"] - valid["fastani_ani"]) * 100
valid["skani_error"] = (valid["skani_ani"] - valid["fastani_ani"]) * 100

for bin_range, sub in valid.groupby("bin"):
    if len(sub) < 5:
        continue
    print(f"\nBin {bin_range}: n={len(sub)}")
    print(f"  FastANI mean: {sub['ani_pct'].mean():.2f}%")
    print(f"  Syn2bANI raw error: mean={sub['s2b_error'].mean():.3f}%, MAE={sub['s2b_error'].abs().mean():.3f}%")
    print(f"  Syn2bANI v4 error:  mean={sub['s2b_v4_error'].mean():.3f}%, MAE={sub['s2b_v4_error'].abs().mean():.3f}%")
    print(f"  skani error:        mean={sub['skani_error'].mean():.3f}%, MAE={sub['skani_error'].abs().mean():.3f}%")
