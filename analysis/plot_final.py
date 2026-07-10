#!/usr/bin/env python3
"""Generate final comparison plots and updated report."""
import csv, pickle, numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BENCH = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI_benchmark_ecoli")
REAL = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI_benchmark_realistic")

def parse_fasta(path):
    seqs = {}
    cid = None; buf = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if cid: seqs[cid] = ''.join(buf)
                cid = line[1:].split()[0]; buf = []
            elif line: buf.append(line)
        if cid: seqs[cid] = ''.join(buf)
    return seqs

# Load GBRT model
with open(BENCH / "syn2bani_gbrt_debias_model.pkl", 'rb') as f:
    model = pickle.load(f)

def simple_debias(raw_ani, af_q, af_r):
    ani_pct = raw_ani * 100.0
    af_min = min(af_q, af_r)
    correction = 0.02 * (100.0 - ani_pct) * (1.0 - af_min)
    return (ani_pct + correction) / 100.0

# Load training data for plotting
with open(BENCH / "gbrt_training_data.csv") as f:
    reader = csv.DictReader(f)
    rows = list(reader)

for r in rows:
    for k in ['raw_ani', 'af_q', 'af_r', 'shared_tags', 'total_q', 'total_r',
              'containment', 'tag_density_q', 'tag_density_r', 'div', 'enzyme_idx',
              'ground_truth_ani', 'simple_debias', 'gbrt_prediction']:
        r[k] = float(r[k])

# Divergence series only (for clean plot)
div_rows = [r for r in rows if r['div'] > 0 and r['enzyme_idx'] == 0]  # BcgI only
div_rows.sort(key=lambda r: r['div'])

# Remove duplicate div values (keep first enzyme)
seen_div = set()
uniq_div = []
for r in div_rows:
    if r['div'] not in seen_div:
        seen_div.add(r['div'])
        uniq_div.append(r)

# Plot 1: Divergence accuracy with all three methods
fig, ax = plt.subplots(1, 1, figsize=(10, 6))

x = [r['div'] * 100 for r in uniq_div]
y_gt = [r['ground_truth_ani'] * 100 for r in uniq_div]
y_raw = [r['raw_ani'] * 100 for r in uniq_div]
y_simple = [r['simple_debias'] * 100 for r in uniq_div]
y_gbrt = [r['gbrt_prediction'] * 100 for r in uniq_div]

ax.plot(x, y_gt, 'o-', color='black', label='Ground Truth', lw=2.5, ms=10, zorder=10)
ax.plot(x, y_raw, 'v--', color='gray', label='Syn2bANI (raw)', lw=1.5, ms=7, alpha=0.6)
ax.plot(x, y_simple, 's:', color='#E84855', label='Simple Debias', lw=2, ms=8)
ax.plot(x, y_gbrt, '^-', color='#2E86AB', label='GBRT Debias', lw=2.5, ms=10)

ax.set_xlabel('Sequence Divergence (%)', fontsize=14)
ax.set_ylabel('ANI (%)', fontsize=14)
ax.set_title('Syn2bANI Debiasing: Raw vs Simple vs GBRT', fontsize=15, fontweight='bold')
ax.legend(fontsize=12, loc='lower left')
ax.grid(True, alpha=0.3)
ax.set_xlim(-0.1, 5.5)

plt.tight_layout()
plt.savefig(BENCH / 'debiasing_comparison.png', dpi=300, bbox_inches='tight')
plt.close()

print(f"Plot saved: {BENCH / 'debiasing_comparison.png'}")

# Plot 2: Error comparison
fig, ax = plt.subplots(1, 1, figsize=(10, 6))

e_raw = [abs(r['raw_ani'] - r['ground_truth_ani']) * 100 for r in uniq_div]
e_simple = [abs(r['simple_debias'] - r['ground_truth_ani']) * 100 for r in uniq_div]
e_gbrt = [abs(r['gbrt_prediction'] - r['ground_truth_ani']) * 100 for r in uniq_div]

x_pos = range(len(x))
width = 0.25
ax.bar([p - width for p in x_pos], e_raw, width, color='gray', alpha=0.5, label='Raw', edgecolor='black')
ax.bar(x_pos, e_simple, width, color='#E84855', alpha=0.8, label='Simple Debias', edgecolor='black')
ax.bar([p + width for p in x_pos], e_gbrt, width, color='#2E86AB', alpha=0.8, label='GBRT Debias', edgecolor='black')

ax.set_xticks(x_pos)
ax.set_xticklabels([f'{v:.2f}' for v in x], rotation=45, ha='right')
ax.set_xlabel('Sequence Divergence (%)', fontsize=14)
ax.set_ylabel('Absolute ANI Error (%)', fontsize=14)
ax.set_title('Debiasing Error Comparison', fontsize=15, fontweight='bold')
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(BENCH / 'debiasing_error_comparison.png', dpi=300, bbox_inches='tight')
plt.close()

print(f"Plot saved: {BENCH / 'debiasing_error_comparison.png'}")

# Plot 3: Realistic MAG results
with open(REAL / "realistic_mag_results.csv") as f:
    reader = csv.DictReader(f)
    real_rows = list(reader)

for r in real_rows:
    r['ani'] = float(r['ani'])
    r['af_q'] = float(r['af_q'])

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Contamination
ax = axes[0]
contam = [r for r in real_rows if r['scenario'] == 'contamination']
contam.sort(key=lambda r: r['param'])
c_x = [r['param'] * 100 for r in contam]
c_y = [r['ani'] * 100 for r in contam]
c_af = [r['af_q'] * 100 for r in contam]
ax.plot(c_x, c_y, 'o-', color='#2E86AB', label='ANI', lw=2, ms=10)
ax.plot(c_x, c_af, 's--', color='#F18F01', label='AF (Query)', lw=2, ms=8)
ax.set_xlabel('Contamination (%)', fontsize=12)
ax.set_ylabel('Value (%)', fontsize=12)
ax.set_title('Contamination Impact', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Assembly error
ax = axes[1]
err = [r for r in real_rows if r['scenario'] == 'assembly_error']
err.sort(key=lambda r: r['param'])
e_x = [r['param'] * 100 for r in err]
e_y = [r['ani'] * 100 for r in err]
e_af = [r['af_q'] * 100 for r in err]
ax.plot(e_x, e_y, 'o-', color='#2E86AB', label='ANI', lw=2, ms=10)
ax.plot(e_x, e_af, 's--', color='#F18F01', label='AF (Query)', lw=2, ms=8)
ax.set_xlabel('Assembly Error (%)', fontsize=12)
ax.set_ylabel('Value (%)', fontsize=12)
ax.set_title('Assembly Error Impact', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(BENCH / 'realistic_mag_results.png', dpi=300, bbox_inches='tight')
plt.close()

print(f"Plot saved: {BENCH / 'realistic_mag_results.png'}")

print("\nAll plots generated successfully!")
