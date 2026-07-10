#!/usr/bin/env python3
"""
Syn2bANI Head-to-Head Benchmark Pipeline
Uses real E. coli genome (NZ_CP026351.1) as reference.
Generates controlled-divergence queries, fragments for MAG simulation.
Implements Python FastANI for comparison.
"""
import csv
import json
import random
import hashlib
from pathlib import Path
import subprocess
import shutil

random.seed(42)

# Paths
ECOLI_REF = Path("/Users/shihuang/test_data/GCF_002953055.1_ASM295305v1_genomic.fna")
BENCH = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI_benchmark_ecoli")
BENCH.mkdir(exist_ok=True)

SYN2BANI = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI/target/release/syn2bani")

def parse_fasta(path):
    """Parse FASTA, return {id: sequence} dict."""
    seqs = {}
    current_id = None
    current_seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    seqs[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            elif line:
                current_seq.append(line)
        if current_id:
            seqs[current_id] = ''.join(current_seq)
    return seqs

def write_fasta(path, seqs):
    """Write {id: sequence} dict to FASTA."""
    with open(path, 'w') as f:
        for seq_id, seq in seqs.items():
            f.write(f">{seq_id}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

def mutate_sequence(seq, snp_rate):
    """Apply random SNPs at given rate, return mutated sequence."""
    seq_list = list(seq)
    mutations = 0
    for i in range(len(seq_list)):
        if random.random() < snp_rate:
            old = seq_list[i]
            choices = [b for b in 'ACGT' if b != old]
            seq_list[i] = random.choice(choices)
            mutations += 1
    return ''.join(seq_list), mutations

def fragment_sequence(seq, n50_target):
    """Break sequence into fragments with exponential distribution achieving target N50."""
    mean_len = n50_target / 0.693
    fragments = []
    pos = 0
    while pos < len(seq):
        frag_len = max(500, int(random.expovariate(1.0 / mean_len)))
        if pos + frag_len > len(seq):
            frag_len = len(seq) - pos
        fragments.append(seq[pos:pos+frag_len])
        pos += frag_len
    return fragments

# ============================================================================
# Step 1: Prepare reference and query genomes from real E. coli
# ============================================================================
print("=" * 60)
print("STEP 1: Preparing real E. coli benchmark data")
print("=" * 60)

ref_seqs = parse_fasta(ECOLI_REF)
ref_id = list(ref_seqs.keys())[0]
ref_seq = ref_seqs[ref_id]
print(f"Reference: {ref_id}, length: {len(ref_seq):,} bp")

# Write clean reference
write_fasta(BENCH / "reference.fasta", {"reference": ref_seq})

# Generate divergence queries
divergence_levels = [0.0005, 0.001, 0.002, 0.005, 0.010, 0.020, 0.030, 0.050]
queries = {}
for div in divergence_levels:
    q_seq, mutations = mutate_sequence(ref_seq, div)
    gt_ani = (len(q_seq) - mutations) / len(q_seq)
    q_name = f"query_div{div:.4f}"
    write_fasta(BENCH / f"{q_name}.fasta", {q_name: q_seq})
    queries[q_name] = {"div": div, "gt_ani": gt_ani, "mutations": mutations}
    print(f"  {q_name}: {mutations:,} mutations, GT_ANI={gt_ani:.6f}")

# Generate N50 fragmentation series (from 2% divergence query)
q_02_seq = parse_fasta(BENCH / "query_div0.0200.fasta")["query_div0.0200"]
gt_02 = queries["query_div0.0200"]["gt_ani"]

n50_levels = [100_000, 50_000, 20_000, 10_000, 5_000, 2_000, 1_000, 500]
for n50 in n50_levels:
    frags = fragment_sequence(q_02_seq, n50)
    seqs = {f"frag_{i}": f for i, f in enumerate(frags)}
    q_name = f"mag_n50_{n50}"
    write_fasta(BENCH / f"{q_name}.fasta", seqs)
    print(f"  {q_name}: {len(frags)} fragments, N50~{n50}")

# Generate completeness series (from N50=10k)
mag_10k_seqs = parse_fasta(BENCH / "mag_n50_10000.fasta")
contig_names = list(mag_10k_seqs.keys())
for comp in [1.0, 0.8, 0.6, 0.5, 0.3, 0.1]:
    n_keep = max(1, int(len(contig_names) * comp))
    kept = {k: mag_10k_seqs[k] for k in contig_names[:n_keep]}
    q_name = f"mag_comp_{comp:.1f}"
    write_fasta(BENCH / f"{q_name}.fasta", kept)
    print(f"  {q_name}: {n_keep}/{len(contig_names)} contigs ({comp*100:.0f}%)")

# Save metadata
with open(BENCH / "metadata.json", 'w') as f:
    json.dump({"reference": ref_id, "length": len(ref_seq), "queries": queries}, f, indent=2)

print(f"\nAll data written to {BENCH}")

# ============================================================================
# Step 2: Run Syn2bANI on all datasets
# ============================================================================
print("\n" + "=" * 60)
print("STEP 2: Running Syn2bANI")
print("=" * 60)

results = []
ref_path = BENCH / "reference.fasta"

for q_file in sorted(BENCH.glob("*.fasta")):
    if q_file.name == "reference.fasta":
        continue
    cmd = [str(SYN2BANI), "dist", str(q_file), str(ref_path)]
    try:
        output = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        lines = output.stdout.strip().split('\n')
        if len(lines) >= 2:
            # Parse TSV output
            headers = lines[0].split('\t')
            values = lines[1].split('\t')
            row = dict(zip(headers, values))
            row['query_name'] = q_file.stem
            row['ani'] = float(row['ani'])
            row['af_q'] = float(row['af_q'])
            row['af_r'] = float(row['af_r'])
            row['shared_tags'] = int(row['shared_tags'])
            results.append(row)
            print(f"  {q_file.stem}: ANI={row['ani']:.4f}, AF={row['af_q']:.4f}, tags={row['shared_tags']}")
    except Exception as e:
        print(f"  ERROR {q_file.stem}: {e}")

# ============================================================================
# Step 3: Implement Python FastANI (k-mer based ANI)
# ============================================================================
print("\n" + "=" * 60)
print("STEP 3: Running Python FastANI")
print("=" * 60)

def fastani_ani(query_seq, ref_seq, k=16, frag_len=3000, min_frag=700):
    """
    Simplified FastANI algorithm:
    1. Fragment query into ~frag_len pieces
    2. Extract k-mers from each fragment
    3. Find matching k-mers in reference
    4. For matching fragments, count shared k-mers / total k-mers
    5. Average over all matching fragments
    """
    # Build reference k-mer index
    ref_kmers = set()
    for i in range(len(ref_seq) - k + 1):
        kmer = ref_seq[i:i+k]
        if 'N' not in kmer:
            ref_kmers.add(kmer)
    
    # Fragment query
    fragments = []
    pos = 0
    while pos < len(query_seq):
        end = min(pos + frag_len, len(query_seq))
        frag = query_seq[pos:end]
        if len(frag) >= min_frag:
            fragments.append(frag)
        pos = end
    
    if not fragments:
        return 0.0, 0.0
    
    matching_frags = 0
    total_ani = 0.0
    
    for frag in fragments:
        frag_kmers = set()
        for i in range(len(frag) - k + 1):
            kmer = frag[i:i+k]
            if 'N' not in kmer:
                frag_kmers.add(kmer)
        
        if not frag_kmers:
            continue
        
        shared = len(frag_kmers & ref_kmers)
        total = len(frag_kmers)
        
        # FastANI uses a threshold: at least 50 shared k-mers for a fragment to count
        if shared >= 50:
            frag_ani = shared / total
            total_ani += frag_ani
            matching_frags += 1
    
    if matching_frags == 0:
        return 0.0, 0.0
    
    avg_ani = total_ani / matching_frags
    af = matching_frags / len(fragments)
    
    # Convert to Mash-style ANI estimate: ANI ≈ (shared_kmers/total_kmers)^(1/k)
    # FastANI uses a more sophisticated regression, but this is a close approximation
    ani_estimate = avg_ani ** (1.0 / k) * 100
    
    return ani_estimate, af

# Run FastANI on all queries
fastani_results = {}
for q_file in sorted(BENCH.glob("*.fasta")):
    if q_file.name == "reference.fasta":
        continue
    q_seqs = parse_fasta(q_file)
    # Concatenate all contigs for comparison
    q_seq = ''.join(q_seqs.values())
    ani, af = fastani_ani(q_seq, ref_seq)
    fastani_results[q_file.stem] = {"ani": ani, "af": af}
    print(f"  {q_file.stem}: ANI={ani:.4f}, AF={af:.4f}")

# ============================================================================
# Step 4: Compute ground truth and compare
# ============================================================================
print("\n" + "=" * 60)
print("STEP 4: Computing ground truth and comparison")
print("=" * 60)

# Calculate ground truth for each query
for row in results:
    q_name = row['query_name']
    q_path = BENCH / f"{q_name}.fasta"
    q_seqs = parse_fasta(q_path)
    q_seq = ''.join(q_seqs.values())
    
    min_len = min(len(ref_seq), len(q_seq))
    matches = sum(1 for i in range(min_len) if ref_seq[i] == q_seq[i])
    gt_ani = matches / min_len
    
    row['ground_truth_ani'] = gt_ani
    row['syn2bani_error'] = abs(row['ani'] - gt_ani) * 100
    
    if q_name in fastani_results:
        row['fastani_ani'] = fastani_results[q_name]['ani'] / 100.0  # normalize to 0-1
        row['fastani_af'] = fastani_results[q_name]['af']
        row['fastani_error'] = abs(row['fastani_ani'] - gt_ani) * 100
    else:
        row['fastani_ani'] = 0.0
        row['fastani_af'] = 0.0
        row['fastani_error'] = 0.0

# Print comparison table
print(f"\n{'Query':<25} {'GT_ANI':<10} {'Syn2bANI':<10} {'S2bErr%':<10} {'FastANI':<10} {'FANIErr%':<10} {'Shared':<8}")
print("-" * 95)
for row in results:
    print(f"{row['query_name']:<25} {row['ground_truth_ani']:<10.4f} {row['ani']:<10.4f} {row['syn2bani_error']:<10.4f} {row['fastani_ani']:<10.4f} {row['fastani_error']:<10.4f} {row['shared_tags']:<8}")

# Save results
with open(BENCH / "comparison_results.csv", 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['query_name', 'ground_truth_ani', 'syn2bani_ani', 'syn2bani_error',
                     'fastani_ani', 'fastani_error', 'syn2bani_af', 'fastani_af', 'shared_tags'])
    for row in results:
        writer.writerow([row['query_name'], row['ground_truth_ani'], row['ani'],
                         row['syn2bani_error'], row['fastani_ani'], row['fastani_error'],
                         row['af_q'], row['fastani_af'], row['shared_tags']])

print(f"\nResults saved to {BENCH / 'comparison_results.csv'}")
