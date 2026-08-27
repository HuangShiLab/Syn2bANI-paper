#!/usr/bin/env python3
"""Analyze alignment fraction (shared fraction) vs FastANI ANI in GTDB-R207 data."""
import pandas as pd
import numpy as np

df = pd.read_csv("/lustre1/g/aos_shihuang/Syn2bANI-paper/results/matrix_gtdb_r207_100k.tsv", sep="\t")
valid = df.dropna(subset=["fastani_ani"]).copy()
valid["ani_pct"] = valid["fastani_ani"] * 100

# Use skani_align_frac as proxy for alignable fraction
# Also compute Syn2bANI af_min
valid["s2b_af_min"] = valid[["s2b_af_q", "s2b_af_r"]].min(axis=1)

print(f"Total valid pairs (with FastANI): {len(valid)}")
print()

bins = [0, 75, 80, 85, 90, 95, 100]
valid["bin"] = pd.cut(valid["ani_pct"], bins=bins)

summary = valid.groupby("bin").agg(
    n=("fastani_ani", "count"),
    mean_ani=("ani_pct", "mean"),
    median_skani_af=("skani_align_frac", "median"),
    mean_skani_af=("skani_align_frac", "mean"),
    median_s2b_af_min=("s2b_af_min", "median"),
    mean_s2b_af_min=("s2b_af_min", "mean"),
    median_s2b_shared=("s2b_shared_tags", "median"),
).reset_index()
print("=== Alignment fraction by ANI bin ===")
print(summary.to_string(index=False, float_format='%.3f'))

print("\n=== Pairs with skani align frac < 50% ===")
low_af = valid[valid["skani_align_frac"] < 0.5]
print(f"Count: {len(low_af)} / {len(valid)} ({len(low_af)/len(valid)*100:.1f}%)")
if len(low_af) > 0:
    print(f"ANI range: {low_af['ani_pct'].min():.2f}% - {low_af['ani_pct'].max():.2f}%")
    print(f"Mean ANI: {low_af['ani_pct'].mean():.2f}%")

print("\n=== ANI threshold where skani align frac drops below 50% ===")
# Find approximate threshold
sorted_df = valid.sort_values("ani_pct")
for threshold in [95, 90, 85, 80, 75, 70]:
    sub = sorted_df[sorted_df["ani_pct"] <= threshold]
    if len(sub) > 0:
        pct_below_50 = (sub["skani_align_frac"] < 0.5).mean() * 100
        print(f"ANI <= {threshold}%: {pct_below_50:.1f}% of pairs have align frac < 50%")
