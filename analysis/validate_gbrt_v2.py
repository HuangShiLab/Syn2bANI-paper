#!/usr/bin/env python3
"""Validate new GBRT v2 model on completely held-out genomes."""
import random, csv, subprocess, json
from pathlib import Path
import numpy as np

random.seed(42)

SYN2B = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI/target/release/syn2bani")
OLD_SYN2B = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI/target/release/syn2bani")  # same binary, but we can test with/without GBRT

# Load new GBRT model (Python) for comparison
with open("/Users/shihuang/Documents/kimi/workspace/Syn2bANI_gbrt_training_v2/gbrt_v2.pkl", 'rb') as f:
    import pickle
    gbrt_v2 = pickle.load(f)

# Held-out test genomes (from new_genomes, not in training set)
TEST_GENOMES = [
    ("GCA_000477835.1", "/Users/shihuang/Downloads/2bRAD同源性/new_genomes/GCA_000477835.1.fna"),
    ("GCA_000521425.1", "/Users/shihuang/Downloads/2bRAD同源性/new_genomes/GCA_000521425.1.fna"),
    ("GCA_002084605.1", "/Users/shihuang/Downloads/2bRAD同源性/new_genomes/GCA_002084605.1.fna"),
    ("GCA_003131265.1", "/Users/shihuang/Downloads/2bRAD同源性/new_genomes/GCA_003131265.1.fna"),
    ("GCA_003243365.1", "/Users/shihuang/Downloads/2bRAD同源性/new_genomes/GCA_003243365.1.fna"),
    ("GCA_005887575.1", "/Users/shihuang/Downloads/2bRAD同源性/new_genomes/GCA_005887575.1.fna"),
    ("GCA_020254085.1", "/Users/shihuang/Downloads/2bRAD同源性/new_genomes/GCA_020254085.1.fna"),
    ("GCA_946476895.1", "/Users/shihuang/Downloads/2bRAD同源性/new_genomes/GCA_946476895.1.fna"),
    ("GCF_005774075.1", "/Users/shihuang/Downloads/2bRAD同源性/new_genomes/GCF_005774075.1.fna"),
    ("GCF_020254005.1", "/Users/shihuang/Downloads/2bRAD同源性/new_genomes/GCF_020254005.1.fna"),
]

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
            for i in range(0, len(seq), 80): f.write(seq[i:i+80]+'\n')

def mutate_sequence(seq, snp_rate):
    seq_list = list(seq)
    mutations = 0
    for i in range(len(seq_list)):
        if random.random() < snp_rate:
            old = seq_list[i]
            choices = [b for b in 'ACGT' if b != old]
            if choices: seq_list[i] = random.choice(choices); mutations += 1
    return ''.join(seq_list), mutations

OUTPUT = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI_gbrt_validation_v2")
OUTPUT.mkdir(exist_ok=True)

print("="*60)
print("GBRT v2 Validation on Held-Out Genomes")
print("="*60)

results = []

for name, ref_path in TEST_GENOMES:
    ref_seqs = parse_fasta(ref_path)
    if not ref_seqs:
        continue
    ref_seq = list(ref_seqs.values())[0]
    ref_len = len(ref_seq)
    ref_gc = sum(1 for b in ref_seq if b in 'GCgc') / ref_len
    
    print(f"\n{name}: {ref_len:,} bp, GC={ref_gc:.3f}")
    
    ref_out = OUTPUT / f"{name}_ref.fasta"
    write_fasta(ref_out, {name: ref_seq})
    
    for div in [0.001, 0.005, 0.01, 0.02, 0.03, 0.05]:
        q_seq, mutations = mutate_sequence(ref_seq, div)
        gt_ani = (len(q_seq) - mutations) / len(q_seq)
        q_name = f"{name}_div{div:.3f}"
        q_path = OUTPUT / f"{q_name}.fasta"
        write_fasta(q_path, {q_name: q_seq})
        
        cmd = [str(SYN2B), "dist", str(q_path), str(ref_out)]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        raw_ani = 0.0; af_q = 0.0; af_r = 0.0; shared = 0
        for line in result.stdout.strip().split('\n'):
            if line.startswith('/'):
                parts = line.split('\t')
                if len(parts) >= 9:
                    raw_ani = float(parts[2])
                    af_q = float(parts[3])
                    af_r = float(parts[4])
                    shared = int(parts[7])
        
        total_q = int(shared / max(af_q, 0.001))
        total_r = int(shared / max(af_r, 0.001))
        max_tags = max(total_q, total_r, 1)
        containment = shared / max_tags
        div_proxy = 1.0 - raw_ani
        
        # Python GBRT v2 prediction
        features = np.array([[raw_ani, af_q, af_r, shared, containment, div_proxy, ref_gc]])
        gbrt_pred = gbrt_v2.predict(features)[0]
        
        # Simple debias
        simple_pred = raw_ani + 0.02 * (1.0 - raw_ani) * (1.0 - min(af_q, af_r))
        
        error_raw = abs(raw_ani - gt_ani) * 100
        error_simple = abs(simple_pred - gt_ani) * 100
        error_gbrt = abs(gbrt_pred - gt_ani) * 100
        
        print(f"  div={div:.3f}: raw_err={error_raw:.4f}% simple_err={error_simple:.4f}% gbrt_v2_err={error_gbrt:.4f}%")
        
        results.append({
            'name': name, 'div': div, 'gt': gt_ani,
            'raw': raw_ani, 'simple': simple_pred, 'gbrt': gbrt_pred,
            'e_raw': error_raw, 'e_simple': error_simple, 'e_gbrt': error_gbrt,
            'ref_gc': ref_gc, 'ref_len': ref_len
        })

# Save
with open(OUTPUT / 'validation_results.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['name', 'div', 'gt', 'raw', 'simple', 'gbrt', 'e_raw', 'e_simple', 'e_gbrt', 'ref_gc', 'ref_len'])
    for r in results:
        writer.writerow([r['name'], r['div'], r['gt'], r['raw'], r['simple'], r['gbrt'],
                         r['e_raw'], r['e_simple'], r['e_gbrt'], r['ref_gc'], r['ref_len']])

# Summary
print("\n" + "="*60)
print("Summary")
print("="*60)

print(f"\n{'Genome':<20} {'Div':<6} {'RawErr':<10} {'SimpleErr':<12} {'GBRTv2Err':<10}")
print("-" * 58)
for r in results:
    print(f"{r['name']:<20} {r['div']:<6.3f} {r['e_raw']:<10.4f} {r['e_simple']:<12.4f} {r['e_gbrt']:<10.4f}")

avg_raw = np.mean([r['e_raw'] for r in results])
avg_simple = np.mean([r['e_simple'] for r in results])
avg_gbrt = np.mean([r['e_gbrt'] for r in results])
max_raw = np.max([r['e_raw'] for r in results])
max_simple = np.max([r['e_simple'] for r in results])
max_gbrt = np.max([r['e_gbrt'] for r in results])

print(f"\n{'Average':<20} {'':<6} {avg_raw:<10.4f} {avg_simple:<12.4f} {avg_gbrt:<10.4f}")
print(f"{'Max':<20} {'':<6} {max_raw:<10.4f} {max_simple:<12.4f} {max_gbrt:<10.4f}")

print(f"\nResults saved to {OUTPUT / 'validation_results.csv'}")
