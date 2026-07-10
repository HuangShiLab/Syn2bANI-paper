#!/usr/bin/env python3
"""Generate final multi-species validation plots."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BENCH = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI_benchmark_multispecies")

with open(BENCH / "multispecies_results.csv") as f:
    reader = csv.DictReader(f)
    rows = list(reader)

for r in rows:
    for k in ['div', 'gt_ani', 'raw_ani', 'simple_ani', 'gbrt_ani',
              'error_raw', 'error_simple', 'error_gbrt']:
        r[k] = float(r[k])

# Group by species
species = sorted(set(r['species'] for r in rows))
divs = sorted(set(r['div'] for r in rows))

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Left: Error comparison by species
ax = axes[0]
x = range(len(species))
width = 0.25

for i, div in enumerate(divs):
    div_rows = [r for r in rows if abs(r['div'] - div) < 0.0001]
    div_rows.sort(key=lambda r: r['species'])
    
    e_raw = [r['error_raw'] for r in div_rows]
    e_simple = [r['error_simple'] for r in div_rows]
    e_gbrt = [r['error_gbrt'] for r in div_rows]
    
    offset = (i - 0.5) * width
    ax.bar([p + offset for p in x], e_gbrt, width, alpha=0.8, label=f'GBRT ({div*100:.1f}%)', edgecolor='black', color='#2E86AB')

ax.set_xticks(x)
ax.set_xticklabels(species, rotation=45, ha='right')
ax.set_xlabel('Species', fontsize=12)
ax.set_ylabel('Absolute ANI Error (%)', fontsize=12)
ax.set_title('GBRT Debiasing: Cross-Species Validation', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, axis='y')

# Right: Average error comparison
ax = axes[1]
methods = ['Raw', 'Simple', 'GBRT']
raw_avg = sum(r['error_raw'] for r in rows) / len(rows)
simple_avg = sum(r['error_simple'] for r in rows) / len(rows)
gbrt_avg = sum(r['error_gbrt'] for r in rows) / len(rows)

bars = ax.bar(methods, [raw_avg, simple_avg, gbrt_avg], color=['gray', '#E84855', '#2E86AB'], alpha=0.8, edgecolor='black')
ax.set_ylabel('Mean Absolute Error (%)', fontsize=12)
ax.set_title('Average Error Across 5 Species', fontsize=13, fontweight='bold')
ax.grid(True, alpha=0.3, axis='y')

# Add value labels on bars
for bar, val in zip(bars, [raw_avg, simple_avg, gbrt_avg]):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
            f'{val:.4f}%', ha='center', va='bottom', fontsize=11, fontweight='bold')

plt.tight_layout()
plt.savefig(BENCH / 'multispecies_validation.png', dpi=300, bbox_inches='tight')
plt.close()

print(f"Plot saved: {BENCH / 'multispecies_validation.png'}")
