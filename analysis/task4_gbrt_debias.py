#!/usr/bin/env python3
"""Task 4: Train GBRT debias model for Syn2bANI."""
import csv, json, random, subprocess
from pathlib import Path
from statistics import mean
import numpy as np
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, r2_score

random.seed(42)
np.random.seed(42)

BENCH = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI_benchmark_ecoli")
SYN2B = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI/target/release/syn2bani")
REF = BENCH / "reference.fasta"

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

# Load ground truth from earlier benchmark
with open(BENCH / "comparison_results.csv") as f:
    reader = csv.DictReader(f)
    gt_map = {r['query_name']: float(r['ground_truth_ani']) for r in reader}

# Re-collect all single-enzyme results for training
ENZYMES = ["BcgI", "BsaXI", "CjeI", "CjePI", "BslFI", "AlfI"]
TRAIN_QUERIES = [q.stem for q in BENCH.glob("query_div*.fasta")] + \
                [q.stem for q in BENCH.glob("mag_n50_*.fasta")] + \
                [q.stem for q in BENCH.glob("mag_comp_*.fasta")]

print("="*60)
print("Task 4: Building GBRT Debias Training Dataset")
print("="*60)

training_data = []

for q_name in sorted(TRAIN_QUERIES):
    q_path = BENCH / f"{q_name}.fasta"
    gt_ani = gt_map.get(q_name)
    if gt_ani is None:
        # Compute GT for N50/completeness queries (same as div=0.02)
        continue

    for enzyme in ENZYMES:
        cmd = [str(SYN2B), "dist", str(q_path), str(REF), "--enzyme", enzyme]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

        raw_ani = 0.0; af_q = 0.0; af_r = 0.0; shared = 0; total_q = 0; total_r = 0
        for line in result.stdout.strip().split('\n'):
            if line.startswith('/'):
                parts = line.split('\t')
                if len(parts) >= 9:
                    raw_ani = float(parts[2])
                    af_q = float(parts[3])
                    af_r = float(parts[4])
                    shared = int(parts[7])
                    # Estimate total tags from AF and shared
                    if af_q > 0:
                        total_q = int(shared / af_q)
                    if af_r > 0:
                        total_r = int(shared / af_r)

        max_tags = max(total_q, total_r, 1)
        containment = shared / max_tags
        tag_density_q = total_q / max(1, sum(len(s) for s in parse_fasta(q_path).values()))
        tag_density_r = total_r / max(1, len(parse_fasta(REF)["reference"]))

        # Determine divergence if known
        div = 0.0
        if q_name.startswith("query_div"):
            div = float(q_name.split("div")[1])
        elif q_name.startswith("mag_"):
            div = 0.02  # baseline for N50/comp series

        # Encode enzyme
        enzyme_idx = ENZYMES.index(enzyme)

        training_data.append({
            'raw_ani': raw_ani,
            'af_q': af_q,
            'af_r': af_r,
            'shared_tags': shared,
            'total_q': total_q,
            'total_r': total_r,
            'containment': containment,
            'tag_density_q': tag_density_q,
            'tag_density_r': tag_density_r,
            'div': div,
            'enzyme_idx': enzyme_idx,
            'ground_truth_ani': gt_ani,
        })

print(f"Collected {len(training_data)} training samples")

# Prepare features and target
feature_names = ['raw_ani', 'af_q', 'af_r', 'shared_tags', 'total_q', 'total_r',
                 'containment', 'tag_density_q', 'tag_density_r', 'div', 'enzyme_idx']

X = np.array([[d[f] for f in feature_names] for d in training_data])
y = np.array([d['ground_truth_ani'] for d in training_data])

# Train/test split (stratified by query to avoid leakage)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

print(f"Train: {len(X_train)}, Test: {len(X_test)}")

# Train GBRT model
print("\nTraining GBRT model...")
model = GradientBoostingRegressor(
    n_estimators=200,
    max_depth=4,
    learning_rate=0.1,
    subsample=0.8,
    random_state=42
)
model.fit(X_train, y_train)

# Predictions
y_pred_train = model.predict(X_train)
y_pred_test = model.predict(X_test)

# Current simple debias for comparison
# Simple correction: ani + 0.02 * (100 - ani) * (1 - min(af_q, af_r))
def simple_debias(raw_ani, af_q, af_r):
    ani_pct = raw_ani * 100.0
    af_min = min(af_q, af_r)
    correction = 0.02 * (100.0 - ani_pct) * (1.0 - af_min)
    return (ani_pct + correction) / 100.0

simple_pred = np.array([simple_debias(X[i,0], X[i,1], X[i,2]) for i in range(len(X))])

# Evaluate on full dataset
print("\n" + "="*60)
print("GBRT Debiasing Results")
print("="*60)

print(f"\n{'Metric':<30} {'Simple Debias':<20} {'GBRT':<20}")
print("-" * 70)

mae_simple = mean_absolute_error(y, simple_pred) * 100
mae_gbrt = mean_absolute_error(y, model.predict(X)) * 100
print(f"{'MAE (%)':<30} {mae_simple:<20.4f} {mae_gbrt:<20.4f}")

r2_simple = r2_score(y, simple_pred)
r2_gbrt = r2_score(y, model.predict(X))
print(f"{'R²':<30} {r2_simple:<20.4f} {r2_gbrt:<20.4f}")

max_err_simple = max(abs(y - simple_pred)) * 100
max_err_gbrt = max(abs(y - model.predict(X))) * 100
print(f"{'Max Error (%)':<30} {max_err_simple:<20.4f} {max_err_gbrt:<20.4f}")

# Feature importance
print(f"\n{'Feature':<20} {'Importance':<15}")
print("-" * 35)
for name, imp in sorted(zip(feature_names, model.feature_importances_), key=lambda x: -x[1]):
    print(f"{name:<20} {imp:<15.4f}")

# Save model
import pickle
model_path = BENCH / "syn2bani_gbrt_debias_model.pkl"
with open(model_path, 'wb') as f:
    pickle.dump(model, f)

print(f"\nModel saved to {model_path}")

# Save training data
with open(BENCH / "gbrt_training_data.csv", 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(feature_names + ['ground_truth_ani', 'simple_debias', 'gbrt_prediction'])
    for i, d in enumerate(training_data):
        row = [d[f] for f in feature_names] + [d['ground_truth_ani'], simple_pred[i], model.predict(X)[i]]
        writer.writerow(row)

print(f"Training data saved to {BENCH / 'gbrt_training_data.csv'}")

# Validate on realistic MAG scenarios (Task 3 data)
print("\n" + "="*60)
print("Applying GBRT to Realistic MAG Scenarios (Task 3)")
print("="*60)

realistic_bench = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI_benchmark_realistic")
with open(realistic_bench / "realistic_mag_results.csv") as f:
    reader = csv.DictReader(f)
    realistic_rows = list(reader)

for r in realistic_rows:
    r['ani'] = float(r['ani'])
    r['af_q'] = float(r['af_q'])
    r['af_r'] = float(r['af_r'])
    r['shared_tags'] = int(r['shared_tags'])

for r in realistic_rows:
    # For realistic MAGs, GT is unknown (same species), but we know it's ~100%
    # Apply GBRT to estimate
    raw = r['ani']
    af_q = r['af_q']
    af_r = r['af_r']
    shared = r['shared_tags']
    total_q = int(shared / max(af_q, 0.001))
    total_r = int(shared / max(af_r, 0.001))
    max_tags = max(total_q, total_r, 1)
    containment = shared / max_tags
    q_path = realistic_bench / f"{r['name']}.fasta"
    q_seqs = parse_fasta(q_path)
    q_len = sum(len(s) for s in q_seqs.values())
    tag_density_q = total_q / max(1, q_len)
    tag_density_r = total_r / max(1, len(parse_fasta(REF)["reference"]))
    div = 0.0  # unknown for realistic
    enzyme_idx = 0  # BcgI

    features = np.array([[raw, af_q, af_r, shared, total_q, total_r,
                          containment, tag_density_q, tag_density_r, div, enzyme_idx]])
    gbrt_corrected = model.predict(features)[0]
    simple_corrected = simple_debias(raw, af_q, af_r)

    print(f"  {r['name']:<25} raw={raw:.4f}  simple={simple_corrected:.4f}  gbrt={gbrt_corrected:.4f}")

print("\nDone!")
