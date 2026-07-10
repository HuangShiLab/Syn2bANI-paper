#!/usr/bin/env python3
"""Export GBRT model to Rust-compatible JSON decision trees."""
import json, pickle, numpy as np
from pathlib import Path

BENCH = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI_benchmark_ecoli")

with open(BENCH / "syn2bani_gbrt_debias_model.pkl", 'rb') as f:
    model = pickle.load(f)

print(f"Model: {type(model).__name__}")
print(f"N_estimators: {model.n_estimators}")
print(f"Max depth: {model.max_depth}")
print(f"Learning rate: {model.learning_rate}")
print(f"Init value: {model.init_.constant_[0][0]}")

feature_names = ['raw_ani', 'af_q', 'af_r', 'shared_tags', 'total_q', 'total_r',
                 'containment', 'tag_density_q', 'tag_density_r', 'div', 'enzyme_idx']

# Export trees as JSON
trees = []
for i in range(model.n_estimators):
    tree = model.estimators_[i, 0].tree_
    n_nodes = tree.node_count
    
    nodes = []
    for node_id in range(n_nodes):
        if tree.children_left[node_id] == tree.children_right[node_id]:
            # Leaf node
            nodes.append({
                'type': 'leaf',
                'value': float(tree.value[node_id][0][0])
            })
        else:
            # Split node
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

json_path = BENCH / "gbrt_model.json"
with open(json_path, 'w') as f:
    json.dump(export, f, indent=2)

print(f"\nModel exported to {json_path}")
print(f"File size: {json_path.stat().st_size / 1024:.1f} KB")

# Verify by running a prediction through the JSON model
def predict_json(features, export_data):
    """Pure Python inference to verify JSON export."""
    prediction = export_data['meta']['init_value']
    lr = export_data['meta']['learning_rate']
    
    for tree in export_data['trees']:
        node_id = 0
        while True:
            node = tree['nodes'][node_id]
            if node['type'] == 'leaf':
                prediction += lr * node['value']
                break
            
            feature_idx = node['feature']
            threshold = node['threshold']
            if features[feature_idx] <= threshold:
                node_id = node['left']
            else:
                node_id = node['right']
    
    return prediction

# Test a few samples
with open(BENCH / "gbrt_training_data.csv") as f:
    import csv
    reader = csv.DictReader(f)
    test_rows = list(reader)[:5]

for r in test_rows:
    features = [float(r[f]) for f in feature_names]
    expected = float(r['gbrt_prediction'])
    json_pred = predict_json(features, export)
    print(f"  Expected: {expected:.6f}, JSON: {json_pred:.6f}, diff: {abs(expected-json_pred):.10f}")

# Also create a compact binary version (for Rust)
# Rust will parse JSON, but we can also create a compact binary format

# Let's also create a simplified model: only top 4 features (div, af_q, raw_ani, af_r)
# Retrain with only these features
from sklearn.ensemble import GradientBoostingRegressor
import csv

with open(BENCH / "gbrt_training_data.csv") as f:
    reader = csv.DictReader(f)
    rows = list(reader)

top_features = ['raw_ani', 'af_q', 'af_r', 'div']
X_top = np.array([[float(r[f]) for f in top_features] for r in rows])
y_top = np.array([float(r['ground_truth_ani']) for r in rows])

model_top = GradientBoostingRegressor(
    n_estimators=200, max_depth=4, learning_rate=0.1, subsample=0.8, random_state=42
)
model_top.fit(X_top, y_top)

y_pred_top = model_top.predict(X_top)
mae_top = np.mean(np.abs(y_pred_top - y_top)) * 100
print(f"\nTop-4 features model MAE: {mae_top:.4f}%")

# Export top-4 model
trees_top = []
for i in range(model_top.n_estimators):
    tree = model_top.estimators_[i, 0].tree_
    n_nodes = tree.node_count
    nodes = []
    for node_id in range(n_nodes):
        if tree.children_left[node_id] == tree.children_right[node_id]:
            nodes.append({'type': 'leaf', 'value': float(tree.value[node_id][0][0])})
        else:
            nodes.append({
                'type': 'split',
                'feature': int(tree.feature[node_id]),
                'feature_name': top_features[int(tree.feature[node_id])],
                'threshold': float(tree.threshold[node_id]),
                'left': int(tree.children_left[node_id]),
                'right': int(tree.children_right[node_id])
            })
    trees_top.append({'nodes': nodes})

export_top = {
    'meta': {
        'n_estimators': model_top.n_estimators,
        'max_depth': model_top.max_depth,
        'learning_rate': model_top.learning_rate,
        'init_value': float(model_top.init_.constant_[0][0]),
        'feature_names': top_features,
    },
    'trees': trees_top
}

json_path_top = BENCH / "gbrt_model_top4.json"
with open(json_path_top, 'w') as f:
    json.dump(export_top, f, indent=2)

print(f"Top-4 model exported to {json_path_top}")
print(f"File size: {json_path_top.stat().st_size / 1024:.1f} KB")

# Also verify top-4 model
for r in test_rows:
    features = [float(r[f]) for f in top_features]
    expected = float(r['gbrt_prediction'])
    json_pred = predict_json(features, export_top)
    print(f"  Top-4 Expected: {expected:.6f}, JSON: {json_pred:.6f}, diff: {abs(expected-json_pred):.10f}")
