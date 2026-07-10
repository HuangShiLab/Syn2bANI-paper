#!/usr/bin/env python3
"""
Simulate structural variations in E. coli genome and validate Syn2bANI detection.
"""
import random
from pathlib import Path
import csv
import subprocess

random.seed(42)

BENCH = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI_benchmark_ecoli")
SYN2BANI = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI/target/release/syn2bani")

def parse_fasta(path):
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
    with open(path, 'w') as f:
        for seq_id, seq in seqs.items():
            f.write(f">{seq_id}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

ref_seqs = parse_fasta(BENCH / "reference.fasta")
ref_seq = ref_seqs["reference"]

print("=" * 60)
print("STEP 1: Simulating Structural Variations")
print("=" * 60)

# Test 1: Inversion (reverse a 50kb segment)
def simulate_inversion(seq, start, end):
    """Reverse a segment of the sequence."""
    before = seq[:start]
    inv = seq[start:end][::-1]
    after = seq[end:]
    return before + inv + after, (start, end)

# Test 2: Translocation (move a 20kb segment to a new position)
def simulate_translocation(seq, start, end, insert_pos):
    """Cut a segment and insert it at a new position."""
    segment = seq[start:end]
    remaining = seq[:start] + seq[end:]
    # Adjust insert position if after the cut region
    if insert_pos > end:
        insert_pos -= (end - start)
    new_seq = remaining[:insert_pos] + segment + remaining[insert_pos:]
    return new_seq, (start, end, insert_pos)

# Test 3: Deletion (remove a 10kb segment)
def simulate_deletion(seq, start, end):
    """Delete a segment."""
    return seq[:start] + seq[end:], (start, end)

# Test 4: Insertion (duplicate a 5kb segment)
def simulate_insertion(seq, start, end, insert_pos):
    """Duplicate a segment and insert elsewhere."""
    segment = seq[start:end]
    return seq[:insert_pos] + segment + seq[insert_pos:], (start, end, insert_pos)

# Test 5: Combined (inversion + deletion + insertion)
def simulate_combined(seq):
    """Apply multiple SVs."""
    sv_log = []
    
    # Inversion at 1Mb
    seq, inv_info = simulate_inversion(seq, 1_000_000, 1_050_000)
    sv_log.append(("inversion", inv_info[0], inv_info[1]))
    
    # Translocation: move 20kb from 2Mb to 3Mb
    seq, trans_info = simulate_translocation(seq, 2_000_000, 2_020_000, 3_000_000)
    sv_log.append(("translocation", trans_info[0], trans_info[1], trans_info[2]))
    
    # Deletion at 500kb
    seq, del_info = simulate_deletion(seq, 500_000, 510_000)
    sv_log.append(("deletion", del_info[0], del_info[1]))
    
    return seq, sv_log

# Generate SV test genomes
sv_tests = {
    "inversion_only": (simulate_inversion(ref_seq, 1_000_000, 1_050_000)[0], [("inversion", 1_000_000, 1_050_000)]),
    "translocation_only": (simulate_translocation(ref_seq, 2_000_000, 2_020_000, 3_000_000)[0], [("translocation", 2_000_000, 2_020_000, 3_000_000)]),
    "deletion_only": (simulate_deletion(ref_seq, 500_000, 510_000)[0], [("deletion", 500_000, 510_000)]),
    "insertion_only": (simulate_insertion(ref_seq, 100_000, 105_000, 200_000)[0], [("insertion", 100_000, 105_000, 200_000)]),
}

# Combined
combined_seq, combined_sv = simulate_combined(ref_seq)
sv_tests["combined"] = (combined_seq, combined_sv)

# Also add SNPs on top of combined (realistic scenario)
snipped = list(combined_seq)
for i in range(len(snipped)):
    if random.random() < 0.002:  # 0.2% SNPs
        old = snipped[i]
        snipped[i] = random.choice([b for b in 'ACGT' if b != old])
sv_tests["combined_with_snps"] = (''.join(snipped), combined_sv + [("snps", 0, 0)])

# Write all test genomes
for name, (seq, sv_log) in sv_tests.items():
    write_fasta(BENCH / f"sv_{name}.fasta", {name: seq})
    print(f"  Generated sv_{name}.fasta: {len(sv_log)} SV events")
    for sv in sv_log:
        print(f"    {sv}")

print("\n" + "=" * 60)
print("STEP 2: Running Syn2bANI with structural analysis")
print("=" * 60)

ref_path = BENCH / "reference.fasta"
sv_results = []

for name, (seq, sv_log) in sv_tests.items():
    q_path = BENCH / f"sv_{name}.fasta"
    
    # Run syn2bani dist for ANI
    cmd_dist = [str(SYN2BANI), "dist", str(q_path), str(ref_path)]
    result_dist = subprocess.run(cmd_dist, capture_output=True, text=True, timeout=60)
    
    # Run syn2bani struct for SV detection
    cmd_struct = [str(SYN2BANI), "struct", str(q_path), str(ref_path), "--rearrangement", "--indel"]
    result_struct = subprocess.run(cmd_struct, capture_output=True, text=True, timeout=60)
    
    # Parse dist output
    dist_ani = 0.0
    dist_af = 0.0
    dist_tags = 0
    for line in result_dist.stdout.strip().split('\n'):
        if line.startswith('/'):
            parts = line.split('\t')
            if len(parts) >= 9:
                dist_ani = float(parts[2])
                dist_af = float(parts[3])
                dist_tags = int(parts[7])
    
    # Parse struct output
    struct_output = result_struct.stdout.strip()
    n_svs_detected = len([l for l in struct_output.split('\n') if l and not l.startswith('query')])
    
    sv_results.append({
        'name': name,
        'true_svs': len(sv_log),
        'detected_svs': n_svs_detected,
        'ani': dist_ani,
        'af': dist_af,
        'shared_tags': dist_tags,
        'struct_output': struct_output,
    })
    
    print(f"\n  {name}:")
    print(f"    True SVs: {len(sv_log)}")
    print(f"    Detected SVs: {n_svs_detected}")
    print(f"    ANI: {dist_ani:.4f}, AF: {dist_af:.4f}")
    if struct_output:
        print(f"    Struct output preview:")
        for line in struct_output.split('\n')[:5]:
            print(f"      {line}")

# Save results
with open(BENCH / "sv_validation_results.csv", 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['name', 'true_svs', 'detected_svs', 'ani', 'af', 'shared_tags'])
    for r in sv_results:
        writer.writerow([r['name'], r['true_svs'], r['detected_svs'], r['ani'], r['af'], r['shared_tags']])

print(f"\nResults saved to {BENCH / 'sv_validation_results.csv'}")
