#!/usr/bin/env python3
"""Generate corrected performance plots from benchmark data."""
import csv, os
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BENCH = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI_performance_benchmark")

with open(BENCH / 'performance_results_corrected.csv') as f:
    rows = list(csv.DictReader(f))

for r in rows:
    for k in r:
        r[k] = float(r[k]) if k != 'n' else int(r[k])

ns = [r['n'] for r in rows]

fig, axes = plt.subplots(2, 3, figsize=(16, 9))

# Row 1: Time
# Sketch time
ax = axes[0, 0]
ax.plot(ns, [r['s2b_sketch_time'] for r in rows], 'o-', color='#2E86AB', label='Syn2bANI', lw=2.5, ms=10)
ax.plot(ns, [r['skani_sketch_time'] for r in rows], 's--', color='#E84855', label='skani', lw=2.5, ms=10)
ax.set_xlabel('Number of Genomes', fontsize=12)
ax.set_ylabel('Sketch Time (s)', fontsize=12)
ax.set_title('A. Sketch Time', fontsize=13, fontweight='bold')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(fontsize=10, loc='upper left')
ax.grid(True, alpha=0.3)

# Query time per pair
ax = axes[0, 1]
s2b_per_pair = [r['s2b_query_time'] / max(1, r['n_pairs']) for r in rows]
skani_per_pair = [r['skani_query_time'] / max(1, r['n_pairs']) for r in rows]
ax.plot(ns, s2b_per_pair, 'o-', color='#2E86AB', label='Syn2bANI', lw=2.5, ms=10)
ax.plot(ns, skani_per_pair, 's--', color='#E84855', label='skani', lw=2.5, ms=10)
ax.set_xlabel('Number of Genomes', fontsize=12)
ax.set_ylabel('Query Time per Pair (s)', fontsize=12)
ax.set_title('B. Query Time (per pair)', fontsize=13, fontweight='bold')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(fontsize=10, loc='upper left')
ax.grid(True, alpha=0.3)

# Total time (sketch + query)
ax = axes[0, 2]
ax.plot(ns, [r['s2b_total_time'] for r in rows], 'o-', color='#2E86AB', label='Syn2bANI', lw=2.5, ms=10)
ax.plot(ns, [r['skani_total_time'] for r in rows], 's--', color='#E84855', label='skani', lw=2.5, ms=10)
ax.set_xlabel('Number of Genomes', fontsize=12)
ax.set_ylabel('Total Time (s)', fontsize=12)
ax.set_title('C. Total Time (sketch + query)', fontsize=13, fontweight='bold')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(fontsize=10, loc='upper left')
ax.grid(True, alpha=0.3)

# Row 2: Memory and Sketch size
# Sketch size per genome
ax = axes[1, 0]
s2b_per_gen = [r['s2b_sketch_size_kb'] / r['n'] for r in rows]
skani_per_gen = [r['skani_sketch_size_kb'] / r['n'] for r in rows]
ax.plot(ns, s2b_per_gen, 'o-', color='#2E86AB', label='Syn2bANI', lw=2.5, ms=10)
ax.plot(ns, skani_per_gen, 's--', color='#E84855', label='skani', lw=2.5, ms=10)
ax.set_xlabel('Number of Genomes', fontsize=12)
ax.set_ylabel('Sketch Size per Genome (KB)', fontsize=12)
ax.set_title('D. Sketch Size (per genome)', fontsize=13, fontweight='bold')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(fontsize=10, loc='upper left')
ax.grid(True, alpha=0.3)

# Total sketch size
ax = axes[1, 1]
ax.plot(ns, [r['s2b_sketch_size_kb'] / 1024 for r in rows], 'o-', color='#2E86AB', label='Syn2bANI', lw=2.5, ms=10)
ax.plot(ns, [r['skani_sketch_size_kb'] / 1024 for r in rows], 's--', color='#E84855', label='skani', lw=2.5, ms=10)
ax.set_xlabel('Number of Genomes', fontsize=12)
ax.set_ylabel('Total Sketch Size (MB)', fontsize=12)
ax.set_title('E. Total Sketch Size', fontsize=13, fontweight='bold')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(fontsize=10, loc='upper left')
ax.grid(True, alpha=0.3)

# Speedup ratio
ax = axes[1, 2]
speedups = [r['skani_total_time'] / r['s2b_total_time'] for r in rows]
ax.plot(ns, speedups, 'D-', color='#2E86AB', lw=2.5, ms=10)
ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5, label='Equal speed')
ax.fill_between(ns, 0, 1, alpha=0.1, color='green', label='Syn2bANI faster')
ax.fill_between(ns, 1, max(speedups)*1.1, alpha=0.1, color='red', label='skani faster')
ax.set_xlabel('Number of Genomes', fontsize=12)
ax.set_ylabel('skani / Syn2bANI Time Ratio', fontsize=12)
ax.set_title('F. Speedup Ratio', fontsize=13, fontweight='bold')
ax.set_xscale('log')
ax.legend(fontsize=10, loc='upper right')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(BENCH / 'performance_comparison_corrected.png', dpi=300, bbox_inches='tight')
plt.close()

print(f"Plot saved: {BENCH / 'performance_comparison_corrected.png'}")

# Also generate a single summary plot (skani Figure 2 style)
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
ax.plot(ns, [r['s2b_total_time'] for r in rows], 'o-', color='#2E86AB', label='Syn2bANI', lw=2.5, ms=12)
ax.plot(ns, [r['skani_total_time'] for r in rows], 's--', color='#E84855', label='skani', lw=2.5, ms=12)
ax.set_xlabel('Number of Genomes', fontsize=14)
ax.set_ylabel('Total Time (s)', fontsize=14)
ax.set_title('Performance Scaling: Syn2bANI vs skani', fontsize=15, fontweight='bold')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(fontsize=12, loc='upper left')
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(BENCH / 'skani_figure2_style.png', dpi=300, bbox_inches='tight')
plt.close()

print(f"Figure 2 style plot saved: {BENCH / 'skani_figure2_style.png'}")

# Print summary
print(f"\n{'N':<5} {'S2B Total':<12} {'Skani Total':<14} {'Speedup':<10} {'S2B KB/gen':<12} {'Skani KB/gen':<14}")
print("-" * 72)
for r in rows:
    speedup = r['skani_total_time'] / r['s2b_total_time'] if r['s2b_total_time'] > 0 else 0
    print(f"{r['n']:<5} {r['s2b_total_time']:<12.2f} {r['skani_total_time']:<14.2f} {speedup:<10.2f} {r['s2b_sketch_size_kb']/r['n']:<12.1f} {r['skani_sketch_size_kb']/r['n']:<14.1f}")
