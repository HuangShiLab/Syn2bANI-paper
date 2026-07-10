#!/usr/bin/env python3
"""Retrain GBRT without 'div' (unknown at runtime), use proxy instead."""
import csv, pickle, numpy as np
from pathlib import Path
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_absolute_error, r2_score
import json

BENCH = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI_benchmark_ecoli")

with open(BENCH / "gbrt_training_data.csv") as f:
    reader = csv.DictReader(f)
    rows = list(reader)

for r in rows:
    for k in r:
        r[k] = float(r[k])

# Proxy for div: 1.0 - raw_ani (approximate)
for r in rows:
    r['div_proxy'] = 1.0 - r['raw_ani']

# Train model WITHOUT div, using proxy instead
# Features available at runtime: raw_ani, af_q, af_r, shared_tags, total_q, total_r, containment, div_proxy
feature_names = ['raw_ani', 'af_q', 'af_r', 'shared_tags', 'containment', 'div_proxy']

X = np.array([[r[f] for f in feature_names] for r in rows])
y = np.array([r['ground_truth_ani'] for r in rows])

model = GradientBoostingRegressor(
    n_estimators=200, max_depth=4, learning_rate=0.1, subsample=0.8, random_state=42
)
model.fit(X, y)

y_pred = model.predict(X)
mae = mean_absolute_error(y, y_pred) * 100
r2 = r2_score(y, y_pred)
print(f"No-div model: MAE={mae:.4f}%, R2={r2:.6f}")

# Export to JSON
trees = []
for i in range(model.n_estimators):
    tree = model.estimators_[i, 0].tree_
    nodes = []
    for node_id in range(tree.node_count):
        if tree.children_left[node_id] == tree.children_right[node_id]:
            nodes.append({'type': 'leaf', 'value': float(tree.value[node_id][0][0])})
        else:
            nodes.append({
                'type': 'split',
                'feature': int(tree.feature[node_id]),
                'feature_name': feature_names[int(tree.feature[node_id])],
                'threshold': float(tree.threshold[node_id]),
                'left': int(tree.children_left[node_id]),
                'right': int(tree.children_right[node_id])
            })
    trees.append({'nodes': nodes})

export = {
    'meta': {
        'n_estimators': model.n_estimators,
        'max_depth': model.max_depth,
        'learning_rate': model.learning_rate,
        'init_value': float(model.init_.constant_[0][0]),
        'feature_names': feature_names,
    },
    'trees': trees
}

json_path = BENCH / "gbrt_model_runtime.json"
with open(json_path, 'w') as f:
    json.dump(export, f, indent=2)

print(f"Runtime model exported to {json_path}")
print(f"File size: {json_path.stat().st_size / 1024:.1f} KB")

# Verify
for r in rows[:5]:
    features = [r[f] for f in feature_names]
    expected = r['gbrt_prediction']
    
    # Reconstruct prediction
    pred = export['meta']['init_value']
    lr = export['meta']['learning_rate']
    for tree in export['trees']:
        node_id = 0
        while True:
            node = tree['nodes'][node_id]
            if node['type'] == 'leaf':
                pred += lr * node['value']
                break
            feat = node['feature']
            if features[feat] <= node['threshold']:
                node_id = node['left']
            else:
                node_id = node['right']
    
    print(f"  Expected: {expected:.6f}, JSON: {pred:.6f}, diff: {abs(expected-pred):.10f}")

# Copy to project
import shutil
shutil.copy(json_path, Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI/gbrt_model_runtime.json"))
print("Copied to project root.")
