#!/usr/bin/env python3
"""
Multi-enzyme consensus ANI benchmark for Syn2bANI.
Tests whether combining multiple Type IIB enzymes improves ANI accuracy.
"""
import random
import subprocess
import csv
from pathlib import Path
from statistics import median

random.seed(42)

BENCH = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI_benchmark_ecoli")
SYN2BANI = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI/target/release/syn2bani")
REF = BENCH / "reference.fasta"

# Enzymes to test
ENZYMES = ["BcgI", "BsaXI", "CjeI", "CjePI", "BslFI", "AlfI"]

# Queries to test (focusing on higher divergence where single-enzyme error is largest)
TEST_QUERIES = [
    "query_div0.0100",  # 1% div
    "query_div0.0200",  # 2% div
    "query_div0.0300",  # 3% div
    "query_div0.0500",  # 5% div
]

print("=" * 70)
print("Multi-Enzyme Consensus ANI Benchmark")
print("=" * 70)

results = []

for q_name in TEST_QUERIES:
    q_path = BENCH / f"{q_name}.fasta"
    
    # Ground truth ANI from previous benchmark
    gt_map = {
        "query_div0.0100": 0.990006,
        "query_div0.0200": 0.980029,
        "query_div0.0300": 0.970024,
        "query_div0.0500": 0.950163,
    }
    gt_ani = gt_map[q_name]
    
    print(f"\n{q_name} (GT_ANI={gt_ani:.4f}):")
    
    # Run each enzyme independently
    enzyme_anis = []
    enzyme_tags = []
    
    for enzyme in ENZYMES:
        cmd = [str(SYN2BANI), "dist", str(q_path), str(REF), "--enzyme", enzyme]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        ani = 0.0
        shared = 0
        for line in result.stdout.strip().split('\n'):
            if line.startswith('/'):
                parts = line.split('\t')
                if len(parts) >= 9:
                    ani = float(parts[2])
                    shared = int(parts[7])
        
        enzyme_anis.append(ani)
        enzyme_tags.append(shared)
        print(f"  {enzyme:10s}: ANI={ani:.4f}, tags={shared}")
    
    # Consensus methods
    # 1. Simple mean
    mean_ani = sum(enzyme_anis) / len(enzyme_anis)
    
    # 2. Weighted mean by shared tag count
    total_tags = sum(enzyme_tags)
    if total_tags > 0:
        weighted_mean = sum(a * t for a, t in zip(enzyme_anis, enzyme_tags)) / total_tags
    else:
        weighted_mean = mean_ani
    
    # 3. Median
    median_ani = median(enzyme_anis)
    
    # 4. Trimmed mean (exclude highest and lowest)
    sorted_anis = sorted(enzyme_anis)
    trimmed_mean = sum(sorted_anis[1:-1]) / max(1, len(sorted_anis) - 2)
    
    print(f"\n  Consensus methods:")
    print(f"    Mean:         {mean_ani:.4f} (error={abs(mean_ani-gt_ani)*100:.3f}%)")
    print(f"    WeightedMean: {weighted_mean:.4f} (error={abs(weighted_mean-gt_ani)*100:.3f}%)")
    print(f"    Median:       {median_ani:.4f} (error={abs(median_ani-gt_ani)*100:.3f}%)")
    print(f"    TrimmedMean:  {trimmed_mean:.4f} (error={abs(trimmed_mean-gt_ani)*100:.3f}%)")
    
    # Best single enzyme
    best_single = min(enzyme_anis, key=lambda a: abs(a - gt_ani))
    best_enzyme = ENZYMES[enzyme_anis.index(best_single)]
    print(f"    Best single ({best_enzyme}): {best_single:.4f} (error={abs(best_single-gt_ani)*100:.3f}%)")
    
    results.append({
        'query': q_name,
        'gt_ani': gt_ani,
        'best_single': best_single,
        'best_enzyme': best_enzyme,
        'mean': mean_ani,
        'weighted_mean': weighted_mean,
        'median': median_ani,
        'trimmed_mean': trimmed_mean,
    })

# Summary
print("\n" + "=" * 70)
print("SUMMARY: Does multi-enzyme consensus improve accuracy?")
print("=" * 70)

print(f"\n{'Query':<20} {'BestSingleErr':<15} {'MeanErr':<15} {'WeightedErr':<15} {'MedianErr':<15} {'TrimmedErr':<15}")
print("-" * 95)
for r in results:
    print(f"{r['query']:<20} {abs(r['best_single']-r['gt_ani'])*100:<15.4f} {abs(r['mean']-r['gt_ani'])*100:<15.4f} {abs(r['weighted_mean']-r['gt_ani'])*100:<15.4f} {abs(r['median']-r['gt_ani'])*100:<15.4f} {abs(r['trimmed_mean']-r['gt_ani'])*100:<15.4f}")

# Save
with open(BENCH / "multienzyme_consensus.csv", 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['query', 'gt_ani', 'best_single', 'best_enzyme', 'mean', 'weighted_mean', 'median', 'trimmed_mean'])
    for r in results:
        writer.writerow([r['query'], r['gt_ani'], r['best_single'], r['best_enzyme'], r['mean'], r['weighted_mean'], r['median'], r['trimmed_mean']])

print(f"\nResults saved to {BENCH / 'multienzyme_consensus.csv'}")
