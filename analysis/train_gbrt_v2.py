#!/usr/bin/env python3
"""Build a large-scale Syn2bANI training dataset from all available genomes."""
import os, random, csv, subprocess, json
from pathlib import Path
from multiprocessing import Pool

random.seed(42)

# Directories to search
SEARCH_DIRS = [
    Path("/Users/shihuang/Downloads/2bRAD同源性/complete_genomes"),
    Path("/Users/shihuang/Downloads/2bRAD同源性/new_genomes"),
    Path("/Users/shihuang/Downloads/genome_260305_2"),
    Path("/Users/shihuang/test_data"),
]

SYN2B = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI/target/release/syn2bani")
OUTPUT = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI_gbrt_training_v2")
OUTPUT.mkdir(exist_ok=True)

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

def fragment_sequence(seq, n50_target):
    mean_len = n50_target / 0.693
    frags = []; pos = 0
    while pos < len(seq):
        fl = max(500, int(random.expovariate(1.0/mean_len)))
        if pos + fl > len(seq): fl = len(seq) - pos
        frags.append(seq[pos:pos+fl]); pos += fl
    return frags

# Step 1: Collect and filter genomes
print("="*60)
print("Step 1: Collecting and filtering genomes")
print("="*60)

all_genomes = {}
for d in SEARCH_DIRS:
    if not d.exists():
        continue
    for f in d.iterdir():
        if f.suffix not in ('.fna', '.fasta', '.fa'):
            continue
        # Skip non-genome files
        fname = f.name.lower()
        if any(x in fname for x in ['sample', 'rep_set', 'seed', 'trim', 'group', 'gut', 'merged', 'reference', 'query', 'contam', 'chim', 'dup', 'err_', 'sv_', 'mag_', 'WecC', 'abfA']):
            continue
        try:
            seqs = parse_fasta(f)
            total_len = sum(len(s) for s in seqs.values())
            n_contigs = len(seqs)
            if total_len >= 1_000_000 and n_contigs <= 5000:  # At least 1 Mb, not too fragmented
                # Concatenate all contigs for a single reference sequence
                combined = ''.join(seqs.values())
                all_genomes[f.stem] = {
                    'path': f,
                    'seq': combined,
                    'length': total_len,
                    'contigs': n_contigs,
                    'gc': sum(1 for b in combined if b in 'GCgc') / total_len
                }
        except Exception as e:
            pass

print(f"Found {len(all_genomes)} usable genomes")
for name, info in sorted(all_genomes.items(), key=lambda x: -x[1]['length']):
    print(f"  {name}: {info['length']:,} bp, {info['contigs']} contigs, GC={info['gc']:.3f}")

# Save genome metadata
with open(OUTPUT / 'genome_metadata.json', 'w') as f:
    json.dump({k: {kk: (str(vv) if isinstance(vv, Path) else (len(vv) if kk == 'seq' else vv)) for kk, vv in v.items()} for k, v in all_genomes.items()}, f, indent=2)

# Step 2: Generate training pairs
print("\n" + "="*60)
print("Step 2: Generating training pairs")
print("="*60)

# For each genome, generate:
# 1. SNP variants at different rates
# 2. Fragmented versions (different N50)
# 3. Low-completeness versions

# Use a subset of genomes for training (to save time, pick ~20 diverse genomes)
genome_names = sorted(all_genomes.keys())
# Select diverse genomes by length (short, medium, long)
by_len = sorted(genome_names, key=lambda n: all_genomes[n]['length'])
# Pick 20 genomes evenly spaced by length
step = max(1, len(by_len) // 20)
selected = [by_len[i*step] for i in range(min(20, len(by_len)))]
# Also ensure E. coli and B. subtilis are included
for special in ['GCF_002953055.1_ASM295305v1_genomic', 'B_subtilis']:
    if special not in selected and special in all_genomes:
        selected.append(special)

selected = list(dict.fromkeys(selected))  # deduplicate, preserve order
print(f"Selected {len(selected)} genomes for training:")
for n in selected:
    info = all_genomes[n]
    print(f"  {n}: {info['length']:,} bp")

# Write selected references
ref_dir = OUTPUT / 'references'
ref_dir.mkdir(exist_ok=True)
for name in selected:
    write_fasta(ref_dir / f'{name}.fasta', {name: all_genomes[name]['seq']})

# Generate variants
variant_dir = OUTPUT / 'variants'
variant_dir.mkdir(exist_ok=True)

div_rates = [0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.05]
n50_levels = [100_000, 20_000, 5_000, 1_000]
comp_levels = [1.0, 0.6, 0.3]

pairs = []  # (ref_name, query_name, query_path, div, n50, comp, gt_ani)

for ref_name in selected:
    ref_seq = all_genomes[ref_name]['seq']
    ref_len = len(ref_seq)

    # Divergence variants
    for div in div_rates:
        q_name = f"{ref_name}_div{div:.4f}"
        q_path = variant_dir / f"{q_name}.fasta"
        if not q_path.exists():
            q_seq, mutations = mutate_sequence(ref_seq, div)
            gt = (len(q_seq) - mutations) / len(q_seq)
            write_fasta(q_path, {q_name: q_seq})
            pairs.append((ref_name, q_name, q_path, div, None, 1.0, gt))
        else:
            pairs.append((ref_name, q_name, q_path, div, None, 1.0, 1.0 - div))  # approximate

    # Fragmentation variants (from 2% divergence base)
    base_seq, _ = mutate_sequence(ref_seq, 0.02)
    for n50 in n50_levels:
        q_name = f"{ref_name}_n50_{n50}"
        q_path = variant_dir / f"{q_name}.fasta"
        if not q_path.exists():
            frags = fragment_sequence(base_seq, n50)
            write_fasta(q_path, {f"frag_{i}": f for i, f in enumerate(frags)})
            gt = 0.98  # approximate
            pairs.append((ref_name, q_name, q_path, 0.02, n50, 1.0, gt))

    # Completeness variants (from N50=5000, 2% div)
    base_frags = fragment_sequence(base_seq, 5000)
    for comp in comp_levels:
        n_keep = max(1, int(len(base_frags) * comp))
        q_name = f"{ref_name}_comp_{comp:.1f}"
        q_path = variant_dir / f"{q_name}.fasta"
        if not q_path.exists():
            kept = {f"frag_{i}": base_frags[i] for i in range(n_keep)}
            write_fasta(q_path, kept)
            gt = 0.98
            pairs.append((ref_name, q_name, q_path, 0.02, 5000, comp, gt))

print(f"Generated {len(pairs)} training pairs")

# Step 3: Run Syn2bANI on all pairs
print("\n" + "="*60)
print("Step 3: Running Syn2bANI on all pairs")
print("="*60)

def run_pair(args):
    ref_name, q_name, q_path, div, n50, comp, gt = args
    ref_path = ref_dir / f"{ref_name}.fasta"
    cmd = [str(SYN2B), "dist", str(q_path), str(ref_path)]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        for line in result.stdout.strip().split('\n'):
            if line.startswith('/'):
                parts = line.split('\t')
                if len(parts) >= 9:
                    raw_ani = float(parts[2])
                    af_q = float(parts[3])
                    af_r = float(parts[4])
                    shared = int(parts[7])
                    return {
                        'ref': ref_name, 'query': q_name, 'div': div, 'n50': n50, 'comp': comp,
                        'gt_ani': gt, 'raw_ani': raw_ani, 'af_q': af_q, 'af_r': af_r,
                        'shared_tags': shared, 'ref_len': all_genomes[ref_name]['length'],
                        'ref_gc': all_genomes[ref_name]['gc']
                    }
    except Exception as e:
        print(f"  ERROR {q_name}: {e}")
    return None

# Run in batches to avoid memory issues
results = []
BATCH = 10
for i in range(0, len(pairs), BATCH):
    batch = pairs[i:i+BATCH]
    print(f"  Batch {i//BATCH + 1}/{(len(pairs)-1)//BATCH + 1}: {i}-{min(i+BATCH, len(pairs))-1}")
    for r in batch:
        res = run_pair(r)
        if res: results.append(res)

print(f"\nCollected {len(results)} results")

# Save training data
with open(OUTPUT / 'training_data_v2.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['ref', 'query', 'div', 'n50', 'comp', 'gt_ani', 'raw_ani', 'af_q', 'af_r', 'shared_tags', 'ref_len', 'ref_gc'])
    for r in results:
        writer.writerow([r['ref'], r['query'], r['div'], r['n50'] or '', r['comp'],
                         r['gt_ani'], r['raw_ani'], r['af_q'], r['af_r'],
                         r['shared_tags'], r['ref_len'], r['ref_gc']])

print(f"Training data saved to {OUTPUT / 'training_data_v2.csv'}")

# Step 4: Train GBRT
print("\n" + "="*60)
print("Step 4: Training GBRT")
print("="*60)

import numpy as np
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, r2_score

for r in results:
    total_q = int(r['shared_tags'] / max(r['af_q'], 0.001))
    total_r = int(r['shared_tags'] / max(r['af_r'], 0.001))
    max_tags = max(total_q, total_r, 1)
    r['total_q'] = total_q
    r['total_r'] = total_r
    r['containment'] = r['shared_tags'] / max_tags
    r['tag_density'] = total_q / max(1, r['ref_len'])
    r['div_proxy'] = 1.0 - r['raw_ani']

feature_cols = ['raw_ani', 'af_q', 'af_r', 'shared_tags', 'containment', 'div_proxy', 'ref_gc']
X = np.array([[r[c] for c in feature_cols] for r in results])
y = np.array([r['gt_ani'] for r in results])

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

model = GradientBoostingRegressor(
    n_estimators=300, max_depth=5, learning_rate=0.05, subsample=0.8, random_state=42
)
model.fit(X_train, y_train)

y_pred = model.predict(X_test)
mae = mean_absolute_error(y_test, y_pred) * 100
r2 = r2_score(y_test, y_pred)
print(f"Test MAE: {mae:.4f}%, R2: {r2:.6f}")

# Full dataset
y_pred_full = model.predict(X)
mae_full = mean_absolute_error(y, y_pred_full) * 100
print(f"Full MAE: {mae_full:.4f}%")

# Export
import pickle, json
with open(OUTPUT / 'gbrt_v2.pkl', 'wb') as f:
    pickle.dump(model, f)

# Export JSON
trees = []
for i in range(model.n_estimators):
    tree = model.estimators_[i, 0].tree_
    nodes = []
    for node_id in range(tree.node_count):
        if tree.children_left[node_id] == tree.children_right[node_id]:
            nodes.append({'type': 'leaf', 'value': float(tree.value[node_id][0][0])})
        else:
            nodes.append({
                'type': 'split', 'feature': int(tree.feature[node_id]),
                'feature_name': feature_cols[int(tree.feature[node_id])],
                'threshold': float(tree.threshold[node_id]),
                'left': int(tree.children_left[node_id]),
                'right': int(tree.children_right[node_id])
            })
    trees.append({'nodes': nodes})

export = {
    'meta': {
        'n_estimators': model.n_estimators, 'max_depth': model.max_depth,
        'learning_rate': model.learning_rate,
        'init_value': float(model.init_.constant_[0][0]),
        'feature_names': feature_cols
    },
    'trees': trees
}
with open(OUTPUT / 'gbrt_v2.json', 'w') as f:
    json.dump(export, f)

print(f"Model exported to {OUTPUT / 'gbrt_v2.json'} ({(OUTPUT / 'gbrt_v2.json').stat().st_size / 1024:.1f} KB)")
print(f"Feature importance:")
for name, imp in sorted(zip(feature_cols, model.feature_importances_), key=lambda x: -x[1]):
    print(f"  {name}: {imp:.4f}")

print("\nDone! Copy gbrt_v2.json to project root and rebuild.")
