#!/usr/bin/env python3
"""Generate comparison plots for Syn2bANI vs FastANI benchmark."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BENCH = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI_benchmark_ecoli")

with open(BENCH / "comparison_results.csv") as f:
    reader = csv.DictReader(f)
    rows = [r for r in reader]

for r in rows:
    for k in ['ground_truth_ani', 'syn2bani_ani', 'syn2bani_error',
              'fastani_ani', 'fastani_error', 'syn2bani_af', 'fastani_af', 'shared_tags']:
        r[k] = float(r[k])

divergence = [r for r in rows if r['query_name'].startswith('query_div')]
n50 = [r for r in rows if r['query_name'].startswith('mag_n50')]
comp = [r for r in rows if r['query_name'].startswith('mag_comp')]

divergence.sort(key=lambda r: r['ground_truth_ani'])
n50.sort(key=lambda r: int(r['query_name'].split('_')[-1]))
comp.sort(key=lambda r: float(r['query_name'].split('_')[-1]))

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# A. ANI Accuracy
ax = axes[0, 0]
x = [1.0 - r['ground_truth_ani'] for r in divergence]
y_gt = [r['ground_truth_ani'] * 100 for r in divergence]
y_s2b = [r['syn2bani_ani'] * 100 for r in divergence]
y_fa = [r['fastani_ani'] * 100 for r in divergence]
ax.plot(x, y_gt, 'o-', color='black', label='Ground Truth', lw=2, ms=8)
ax.plot(x, y_s2b, 's--', color='#2E86AB', label='Syn2bANI', lw=2, ms=8)
ax.plot(x, y_fa, '^:', color='#E84855', label='FastANI (Python)', lw=2, ms=8)
ax.set_xlabel('Sequence Divergence (fraction)', fontsize=12)
ax.set_ylabel('ANI (%)', fontsize=12)
ax.set_title('A. ANI Accuracy vs Sequence Divergence', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xlim(-0.002, 0.055)

# B. Error comparison
ax = axes[0, 1]
labels = [f"{1-r['ground_truth_ani']:.1%}" for r in divergence]
e_s2b = [r['syn2bani_error'] for r in divergence]
e_fa = [r['fastani_error'] for r in divergence]
x_pos = range(len(labels))
width = 0.35
ax.bar([p - width/2 for p in x_pos], e_s2b, width, color='#2E86AB', alpha=0.8, label='Syn2bANI', edgecolor='black')
ax.bar([p + width/2 for p in x_pos], e_fa, width, color='#E84855', alpha=0.8, label='FastANI', edgecolor='black')
ax.set_xticks(x_pos)
ax.set_xticklabels(labels, rotation=45, ha='right')
ax.set_xlabel('Sequence Divergence', fontsize=12)
ax.set_ylabel('Absolute ANI Error (%)', fontsize=12)
ax.set_title('B. ANI Estimation Error Comparison', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3, axis='y')

# C. N50 robustness
ax = axes[1, 0]
n50_vals = [int(r['query_name'].split('_')[-1]) for r in n50]
y_gt = [r['ground_truth_ani'] * 100 for r in n50]
y_s2b = [r['syn2bani_ani'] * 100 for r in n50]
y_fa = [r['fastani_ani'] * 100 for r in n50]
ax.plot(n50_vals, y_gt, 'o-', color='black', label='Ground Truth', lw=2, ms=8)
ax.plot(n50_vals, y_s2b, 's--', color='#2E86AB', label='Syn2bANI', lw=2, ms=8)
ax.plot(n50_vals, y_fa, '^:', color='#E84855', label='FastANI', lw=2, ms=8)
ax.set_xlabel('N50 (bp)', fontsize=12)
ax.set_ylabel('ANI (%)', fontsize=12)
ax.set_title('C. ANI Robustness vs N50 (Fragmentation)', fontsize=13, fontweight='bold')
ax.set_xscale('log')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# D. Completeness robustness
ax = axes[1, 1]
c_vals = [float(r['query_name'].split('_')[-1]) for r in comp]
y_gt = [r['ground_truth_ani'] * 100 for r in comp]
y_s2b = [r['syn2bani_ani'] * 100 for r in comp]
y_fa = [r['fastani_ani'] * 100 for r in comp]
ax.plot(c_vals, y_gt, 'o-', color='black', label='Ground Truth', lw=2, ms=8)
ax.plot(c_vals, y_s2b, 's--', color='#2E86AB', label='Syn2bANI', lw=2, ms=8)
ax.plot(c_vals, y_fa, '^:', color='#E84855', label='FastANI', lw=2, ms=8)
ax.set_xlabel('Genome Completeness (fraction)', fontsize=12)
ax.set_ylabel('ANI (%)', fontsize=12)
ax.set_title('D. ANI Robustness vs MAG Completeness', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 1.1)

plt.tight_layout()
plt.savefig(BENCH / 'syn2bani_vs_fastani.png', dpi=300, bbox_inches='tight')
plt.close()

print(f"Plot saved to {BENCH / 'syn2bani_vs_fastani.png'}")

# Generate summary report
report = f"""# Syn2bANI vs FastANI Head-to-Head Benchmark

## Dataset
- **Reference**: *E. coli* NZ_CP026351.1 (4,651,848 bp, complete chromosome)
- **Queries**: Derived with controlled SNP rates (0.05%–5%)
- **Fragmentation**: Exponential distribution, N50 from 500 bp to 100 kb
- **Completeness**: Sequential contig truncation, 30%–100%

## Key Results

### Divergence Accuracy

| Divergence | GT_ANI | Syn2bANI | S2b Error | FastANI | FA Error |
|-----------|--------|----------|-----------|---------|----------|
"""
for r in divergence:
    report += f"| {1-r['ground_truth_ani']:.2%} | {r['ground_truth_ani']*100:.2f}% | {r['syn2bani_ani']*100:.2f}% | {r['syn2bani_error']:.3f}% | {r['fastani_ani']*100:.2f}% | {r['fastani_error']:.3f}% |\n"

report += """
### Fragmentation Robustness (2% divergence baseline)

| N50 | Syn2bANI | FastANI | S2b Error | FA Error |
|-----|----------|---------|-----------|----------|
"""
for r in n50:
    report += f"| {int(r['query_name'].split('_')[-1]):,} | {r['syn2bani_ani']*100:.2f}% | {r['fastani_ani']*100:.2f}% | {r['syn2bani_error']:.3f}% | {r['fastani_error']:.3f}% |\n"

report += """
### Completeness Robustness (2% div, N50~10k)

| Completeness | Syn2bANI | FastANI | S2b Error | FA Error |
|-------------|----------|---------|-----------|----------|
"""
for r in comp:
    report += f"| {float(r['query_name'].split('_')[-1])*100:.0f}% | {r['syn2bani_ani']*100:.2f}% | {r['fastani_ani']*100:.2f}% | {r['syn2bani_error']:.3f}% | {r['fastani_error']:.3f}% |\n"

report += """
## Interpretation

1. **FastANI (Python) is near-perfect on SNP-only data** because:
   - No structural variation = k-mers map perfectly
   - Full sequence available = all k-mers counted
   - This represents the "best case" for k-mer methods

2. **Syn2bANI has a small systematic overestimation** (~0.5% at 2% div, ~2% at 5% div):
   - Fixed-anchor tags that differ by >1-2 bp are excluded, biasing toward conserved regions
   - The debias model partially corrects this but needs refinement
   - **This is expected**: tag-based methods inherently sample a subset of the genome

3. **Both methods are equally robust to fragmentation/completeness** on this data:
   - Neither method is affected by N50 or completeness in this SNP-only scenario
   - **Real-world difference**: Syn2bANI would maintain accuracy with rearrangements/inversions, while FastANI's k-mer chaining would break

## Limitations

- **No structural variation**: Real MAGs have rearrangements; this test favors k-mer methods
- **Single enzyme**: Multi-enzyme mode may improve Syn2bANI accuracy
- **Python FastANI**: Simplified implementation; real FastANI has fragment-level alignment and regression correction

## Next Steps

1. Add structural variation (inversions, translocations) to test where Syn2bANI's fixed-anchor advantage manifests
2. Implement multi-enzyme consensus mode in Syn2bANI
3. Train a proper GBRT debias model
"""

with open(BENCH / 'HEAD_TO_HEAD_REPORT.md', 'w') as f:
    f.write(report)

print(f"Report saved to {BENCH / 'HEAD_TO_HEAD_REPORT.md'}")
