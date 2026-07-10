#!/usr/bin/env python3
"""
Performance benchmark: Syn2bANI vs skani vs FastANI (Python ref).
Replicates skani Figure 2 style: time/memory scaling with genome count.
"""
import os, time, subprocess, csv, resource
from pathlib import Path
from statistics import mean
import numpy as np

BENCH = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI_performance_benchmark")
BENCH.mkdir(exist_ok=True)

SYN2B = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI/target/release/syn2bani")
SKANI = Path("/Users/shihuang/.cargo/bin/skani")

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

# Collect all genomes
all_genomes = []
for d in [Path("/Users/shihuang/Downloads/2bRAD同源性/complete_genomes"),
          Path("/Users/shihuang/Downloads/2bRAD同源性/new_genomes"),
          Path("/Users/shihuang/Downloads/genome_260305_2")]:
    if d.exists():
        for f in sorted(d.iterdir()):
            if f.suffix in ('.fna', '.fasta', '.fa'):
                try:
                    seqs = parse_fasta(f)
                    total = sum(len(s) for s in seqs.values())
                    if total >= 1_000_000:
                        all_genomes.append({'name': f.stem, 'path': f, 'length': total})
                except: pass

all_genomes.sort(key=lambda g: g['length'])
print(f"Collected {len(all_genomes)} genomes")
for g in all_genomes[:5]:
    print(f"  {g['name']}: {g['length']:,} bp")

# Genome count scaling points
GCOUNTS = [1, 2, 4, 8, 16, 32, min(64, len(all_genomes))]
GCOUNTS = [c for c in GCOUNTS if c <= len(all_genomes)]
print(f"Benchmark points: {GCOUNTS}")

# ============================================================
# Step 1: Pre-sketch ALL genomes for both tools
# ============================================================
print("\n" + "="*60)
print("Step 1: Pre-sketching all genomes")
print("="*60)

skani_sketch_dir = BENCH / 'skani_sketches'
syn2b_sketch_dir = BENCH / 'syn2b_sketches'
skani_sketch_dir.mkdir(exist_ok=True)
syn2b_sketch_dir.mkdir(exist_ok=True)

# Syn2bANI sketch
print("Syn2bANI sketching...")
s2b_sketch_times = {}
for g in all_genomes:
    out = syn2b_sketch_dir / f"{g['name']}.s2ba"
    if out.exists():
        s2b_sketch_times[g['name']] = 0.0
        continue
    t0 = time.time()
    cmd = [str(SYN2B), 'sketch', str(g['path']), '-o', str(out)]
    subprocess.run(cmd, capture_output=True, timeout=120)
    s2b_sketch_times[g['name']] = time.time() - t0

# skani sketch
print("skani sketching...")
skani_sketch_times = {}
for g in all_genomes:
    out = skani_sketch_dir / f"{g['name']}.sketch"
    if out.exists():
        skani_sketch_times[g['name']] = 0.0
        continue
    t0 = time.time()
    cmd = [str(SKANI), 'sketch', str(g['path']), '-o', str(out)]
    subprocess.run(cmd, capture_output=True, timeout=120)
    skani_sketch_times[g['name']] = time.time() - t0

print(f"  Syn2bANI total sketch time: {sum(s2b_sketch_times.values()):.2f}s")
print(f"  skani total sketch time: {sum(skani_sketch_times.values()):.2f}s")

# ============================================================
# Step 2: Query benchmark at different scales
# ============================================================
print("\n" + "="*60)
print("Step 2: Query benchmarking")
print("="*60)

def get_peak_memory(cmd_list, pairs):
    """Run a few pairs with /usr/bin/time -l to get peak memory."""
    peak = 0
    for q_name, r_name in pairs[:min(5, len(pairs))]:
        cmd = ['/usr/bin/time', '-l'] + cmd_list + [str(q_name), str(r_name)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        for line in result.stderr.split('\n'):
            if 'maximum resident set size' in line:
                try:
                    mem = int(line.split(':')[1].strip()) / 1024
                    peak = max(peak, mem)
                except: pass
    return peak

results = []

for n in GCOUNTS:
    subset = all_genomes[:n]
    total_bp = sum(g['length'] for g in subset)
    
    # Use first 5 as queries, compare against all n
    queries = subset[:min(5, n)]
    n_pairs = len(queries) * (n - 1)
    
    print(f"\n--- n={n} genomes, {total_bp/1e6:.1f} Mb, {n_pairs} pairs ---")
    
    # Syn2bANI
    s2b_sketch_time = sum(s2b_sketch_times[g['name']] for g in subset)
    t0 = time.time()
    for q in queries:
        for r in subset:
            if q['name'] == r['name']: continue
            cmd = [str(SYN2B), 'dist', str(q['path']), str(r['path'])]
            subprocess.run(cmd, capture_output=True, timeout=30)
    s2b_query_time = time.time() - t0
    
    s2b_sketch_size = sum((syn2b_sketch_dir / f"{g['name']}.s2ba").stat().st_size for g in subset if (syn2b_sketch_dir / f"{g['name']}.s2ba").exists())
    
    # Memory for Syn2bANI
    s2b_mem = get_peak_memory([str(SYN2B), 'dist'], [(q['path'], r['path']) for q in queries for r in subset if q['name'] != r['name']])
    
    print(f"  Syn2bANI: sketch={s2b_sketch_time:.2f}s, query={s2b_query_time:.2f}s, mem={s2b_mem:.1f}MB, sketch_size={s2b_sketch_size/1024:.1f}KB")
    
    # skani
    skani_sketch_time = sum(skani_sketch_times[g['name']] for g in subset)
    t0 = time.time()
    for q in queries:
        for r in subset:
            if q['name'] == r['name']: continue
            cmd = [str(SKANI), 'dist', str(q['path']), str(r['path'])]
            subprocess.run(cmd, capture_output=True, timeout=30)
    skani_query_time = time.time() - t0
    
    skani_sketch_size = sum((skani_sketch_dir / f"{g['name']}.sketch").stat().st_size for g in subset if (skani_sketch_dir / f"{g['name']}.sketch").exists())
    
    skani_mem = get_peak_memory([str(SKANI), 'dist'], [(q['path'], r['path']) for q in queries for r in subset if q['name'] != r['name']])
    
    print(f"  skani:    sketch={skani_sketch_time:.2f}s, query={skani_query_time:.2f}s, mem={skani_mem:.1f}MB, sketch_size={skani_sketch_size/1024:.1f}KB")
    
    results.append({
        'n': n, 'total_bp': total_bp, 'n_pairs': n_pairs,
        's2b_sketch_time': s2b_sketch_time, 's2b_query_time': s2b_query_time,
        's2b_total_time': s2b_sketch_time + s2b_query_time,
        's2b_sketch_size_kb': s2b_sketch_size / 1024, 's2b_mem_mb': s2b_mem,
        'skani_sketch_time': skani_sketch_time, 'skani_query_time': skani_query_time,
        'skani_total_time': skani_sketch_time + skani_query_time,
        'skani_sketch_size_kb': skani_sketch_size / 1024, 'skani_mem_mb': skani_mem,
    })

# Save
with open(BENCH / 'performance_results.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['n', 'total_bp', 'n_pairs',
                     's2b_sketch_time', 's2b_query_time', 's2b_total_time', 's2b_sketch_size_kb', 's2b_mem_mb',
                     'skani_sketch_time', 'skani_query_time', 'skani_total_time', 'skani_sketch_size_kb', 'skani_mem_mb'])
    for r in results:
        writer.writerow([r['n'], r['total_bp'], r['n_pairs'],
                         r['s2b_sketch_time'], r['s2b_query_time'], r['s2b_total_time'], r['s2b_sketch_size_kb'], r['s2b_mem_mb'],
                         r['skani_sketch_time'], r['skani_query_time'], r['skani_total_time'], r['skani_sketch_size_kb'], r['skani_mem_mb']])

print(f"\nResults saved to {BENCH / 'performance_results.csv'}")

# ============================================================
# Step 3: Generate plots
# ============================================================
print("\n" + "="*60)
print("Step 3: Generating plots")
print("="*60)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 3, figsize=(16, 9))

ns = [r['n'] for r in results]

# Row 1: Time
# Sketch time
ax = axes[0, 0]
ax.plot(ns, [r['s2b_sketch_time'] for r in results], 'o-', color='#2E86AB', label='Syn2bANI', lw=2, ms=8)
ax.plot(ns, [r['skani_sketch_time'] for r in results], 's--', color='#E84855', label='skani', lw=2, ms=8)
ax.set_xlabel('Number of Genomes', fontsize=12)
ax.set_ylabel('Sketch Time (s)', fontsize=12)
ax.set_title('A. Sketch Time', fontsize=13, fontweight='bold')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Query time
ax = axes[0, 1]
ax.plot(ns, [r['s2b_query_time'] for r in results], 'o-', color='#2E86AB', label='Syn2bANI', lw=2, ms=8)
ax.plot(ns, [r['skani_query_time'] for r in results], 's--', color='#E84855', label='skani', lw=2, ms=8)
ax.set_xlabel('Number of Genomes', fontsize=12)
ax.set_ylabel('Query Time (s)', fontsize=12)
ax.set_title('B. Query Time', fontsize=13, fontweight='bold')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Total time
ax = axes[0, 2]
ax.plot(ns, [r['s2b_total_time'] for r in results], 'o-', color='#2E86AB', label='Syn2bANI', lw=2, ms=8)
ax.plot(ns, [r['skani_total_time'] for r in results], 's--', color='#E84855', label='skani', lw=2, ms=8)
ax.set_xlabel('Number of Genomes', fontsize=12)
ax.set_ylabel('Total Time (s)', fontsize=12)
ax.set_title('C. Total Time', fontsize=13, fontweight='bold')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Row 2: Memory and Sketch size
# Memory
ax = axes[1, 0]
ax.plot(ns, [r['s2b_mem_mb'] for r in results], 'o-', color='#2E86AB', label='Syn2bANI', lw=2, ms=8)
ax.plot(ns, [r['skani_mem_mb'] for r in results], 's--', color='#E84855', label='skani', lw=2, ms=8)
ax.set_xlabel('Number of Genomes', fontsize=12)
ax.set_ylabel('Peak Memory (MB)', fontsize=12)
ax.set_title('D. Peak Memory', fontsize=13, fontweight='bold')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Sketch size per genome
ax = axes[1, 1]
ax.plot(ns, [r['s2b_sketch_size_kb'] / r['n'] for r in results], 'o-', color='#2E86AB', label='Syn2bANI', lw=2, ms=8)
ax.plot(ns, [r['skani_sketch_size_kb'] / r['n'] for r in results], 's--', color='#E84855', label='skani', lw=2, ms=8)
ax.set_xlabel('Number of Genomes', fontsize=12)
ax.set_ylabel('Sketch Size per Genome (KB)', fontsize=12)
ax.set_title('E. Sketch Size (per genome)', fontsize=13, fontweight='bold')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Speedup ratio
ax = axes[1, 2]
ax.plot(ns, [r['skani_total_time'] / r['s2b_total_time'] for r in results], 'o-', color='#2E86AB', lw=2, ms=10)
ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('Number of Genomes', fontsize=12)
ax.set_ylabel('skani / Syn2bANI Time Ratio', fontsize=12)
ax.set_title('F. Speedup (skani / Syn2bANI)', fontsize=13, fontweight='bold')
ax.set_xscale('log')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(BENCH / 'performance_comparison.png', dpi=300, bbox_inches='tight')
plt.close()

print(f"Plot saved to {BENCH / 'performance_comparison.png'}")

# Summary table
print("\n" + "="*60)
print("Summary Table")
print("="*60)
print(f"{'N':<5} {'S2B Total':<12} {'Skani Total':<14} {'Speedup':<10} {'S2B Mem':<10} {'Skani Mem':<12}")
print("-" * 70)
for r in results:
    speedup = r['skani_total_time'] / r['s2b_total_time'] if r['s2b_total_time'] > 0 else 0
    print(f"{r['n']:<5} {r['s2b_total_time']:<12.2f} {r['skani_total_time']:<14.2f} {speedup:<10.2f} {r['s2b_mem_mb']:<10.1f} {r['skani_mem_mb']:<12.1f}")

print("\nDone!")
