#!/usr/bin/env python3
"""Multi-species validation: test GBRT on 5 different bacterial genomes."""
import random, csv, subprocess
from pathlib import Path
import pickle
import numpy as np

random.seed(42)

BENCH = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI_benchmark_multispecies")
BENCH.mkdir(exist_ok=True)
SYN2B = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI/target/release/syn2bani")

# Load GBRT model
with open("/Users/shihuang/Documents/kimi/workspace/Syn2bANI_benchmark_ecoli/syn2bani_gbrt_debias_model.pkl", 'rb') as f:
    gbrt_model = pickle.load(f)

feature_names = ['raw_ani', 'af_q', 'af_r', 'shared_tags', 'total_q', 'total_r',
                 'containment', 'tag_density_q', 'tag_density_r', 'div', 'enzyme_idx']

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

def write_fasta(path, seqs):
    with open(path, 'w') as f:
        for sid, seq in seqs.items():
            f.write(f'>{sid}\n')
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

def mutate_sequence(seq, snp_rate):
    seq_list = list(seq)
    mutations = 0
    for i in range(len(seq_list)):
        if random.random() < snp_rate:
            old = seq_list[i]
            seq_list[i] = random.choice([b for b in 'ACGT' if b != old])
            mutations += 1
    return ''.join(seq_list), mutations

# Select 5 genomes
SPECIES = [
    ("E_coli", "/Users/shihuang/test_data/GCF_002953055.1_ASM295305v1_genomic.fna"),
    ("B_subtilis", "/Users/shihuang/Downloads/2bRAD同源性/B_subtilis.fasta"),
    ("BB006", "/Users/shihuang/Downloads/genome_260305_2/BB006.fasta"),
    ("BB18", "/Users/shihuang/Downloads/genome_260305_2/BB18.fasta"),
    ("LA100", "/Users/shihuang/Downloads/genome_260305_2/LA100.fasta"),
]

print("="*60)
print("Multi-Species GBRT Validation")
print("="*60)

results = []

for species_name, ref_path in SPECIES:
    ref_seqs = parse_fasta(ref_path)
    if not ref_seqs:
        print(f"  SKIP {species_name}: empty file")
        continue
    ref_id = list(ref_seqs.keys())[0]
    ref_seq = ref_seqs[ref_id]
    
    print(f"\n{species_name}: {ref_id}, {len(ref_seq):,} bp")
    
    # Write reference
    ref_out = BENCH / f"{species_name}_ref.fasta"
    write_fasta(ref_out, {f"{species_name}_ref": ref_seq})
    
    # Generate 2 divergence variants
    for div in [0.001, 0.020]:
        q_seq, mutations = mutate_sequence(ref_seq, div)
        gt_ani = (len(q_seq) - mutations) / len(q_seq)
        q_name = f"{species_name}_div{div:.4f}"
        q_path = BENCH / f"{q_name}.fasta"
        write_fasta(q_path, {q_name: q_seq})
        
        # Run Syn2bANI
        cmd = [str(SYN2B), "dist", str(q_path), str(ref_out)]
        output = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        raw_ani = 0.0; af_q = 0.0; af_r = 0.0; shared = 0; total_q = 0; total_r = 0
        for line in output.stdout.strip().split('\n'):
            if line.startswith('/'):
                parts = line.split('\t')
                if len(parts) >= 9:
                    raw_ani = float(parts[2])
                    af_q = float(parts[3])
                    af_r = float(parts[4])
                    shared = int(parts[7])
                    if af_q > 0: total_q = int(shared / af_q)
                    if af_r > 0: total_r = int(shared / af_r)
        
        max_tags = max(total_q, total_r, 1)
        containment = shared / max_tags
        tag_density_q = total_q / max(1, len(q_seq))
        tag_density_r = total_r / max(1, len(ref_seq))
        
        # GBRT prediction (full model with div)
        features = np.array([[raw_ani, af_q, af_r, shared, total_q, total_r,
                              containment, tag_density_q, tag_density_r, div, 0]])
        gbrt_pred = gbrt_model.predict(features)[0]
        
        # Simple debias
        simple_pred = raw_ani + 0.02 * (1.0 - raw_ani) * (1.0 - min(af_q, af_r))
        
        error_raw = abs(raw_ani - gt_ani) * 100
        error_simple = abs(simple_pred - gt_ani) * 100
        error_gbrt = abs(gbrt_pred - gt_ani) * 100
        
        print(f"  {q_name}: raw={raw_ani:.4f} simple={simple_pred:.4f} gbrt={gbrt_pred:.4f} gt={gt_ani:.4f}")
        print(f"    error: raw={error_raw:.4f}% simple={error_simple:.4f}% gbrt={error_gbrt:.4f}%")
        
        results.append({
            'species': species_name,
            'div': div,
            'gt_ani': gt_ani,
            'raw_ani': raw_ani,
            'simple_ani': simple_pred,
            'gbrt_ani': gbrt_pred,
            'error_raw': error_raw,
            'error_simple': error_simple,
            'error_gbrt': error_gbrt,
            'shared_tags': shared,
        })

# Save
with open(BENCH / "multispecies_results.csv", 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['species', 'div', 'gt_ani', 'raw_ani', 'simple_ani', 'gbrt_ani',
                     'error_raw', 'error_simple', 'error_gbrt', 'shared_tags'])
    for r in results:
        writer.writerow([r['species'], r['div'], r['gt_ani'], r['raw_ani'], r['simple_ani'],
                         r['gbrt_ani'], r['error_raw'], r['error_simple'], r['error_gbrt'],
                         r['shared_tags']])

# Summary
print("\n" + "="*60)
print("Summary")
print("="*60)

print(f"\n{'Species':<15} {'Div':<8} {'RawErr':<10} {'SimpleErr':<12} {'GBRTErr':<10}")
print("-" * 55)
for r in results:
    print(f"{r['species']:<15} {r['div']:<8.3f} {r['error_raw']:<10.4f} {r['error_simple']:<12.4f} {r['error_gbrt']:<10.4f}")

avg_raw = np.mean([r['error_raw'] for r in results])
avg_simple = np.mean([r['error_simple'] for r in results])
avg_gbrt = np.mean([r['error_gbrt'] for r in results])
print(f"\n{'Average':<15} {'':<8} {avg_raw:<10.4f} {avg_simple:<12.4f} {avg_gbrt:<10.4f}")

print(f"\nResults saved to {BENCH / 'multispecies_results.csv'}")
