# Syn2bANI GBRT Debiasing: Technical Deep Dive

> A detailed explanation of why Syn2bANI needs debiasing, how the Gradient Boosted Regression Tree (GBRT) model corrects systematic ANI overestimation, and how it is embedded into the Rust binary.

---

## Table of Contents

1. [The Problem: Systematic ANI Overestimation in Fixed-Anchor Tag Matching](#1-the-problem)
2. [Why Traditional Debiasing Fails](#2-traditional-debiasing)
3. [GBRT: Architecture and Training](#3-gbrt-architecture)
4. [Feature Engineering](#4-feature-engineering)
5. [Model Export: From Python to Rust](#5-model-export)
6. [Rust Inference Engine](#6-rust-inference)
7. [Validation Results](#7-validation)
8. [Usage in Syn2bANI](#8-usage)

---

## 1. The Problem: Systematic ANI Overestimation in Fixed-Anchor Tag Matching

### 1.1 How Syn2bANI Computes ANI

Syn2bANI extracts **2bRAD tags** (~32 bp DNA fragments flanking Type IIB restriction sites) from both query and reference genomes. These tags act as **fixed positional anchors** — unlike random k-mers, each tag has a deterministic location defined by the enzyme recognition sequence.

The raw ANI is computed as:

```
raw_ANI = mean(local_ANI of all matched tag pairs)
```

where `local_ANI` for a single matched pair is:

```
local_ANI = 1 - (Hamming_distance / 32)
```

### 1.2 The Systematic Bias

This approach suffers from a **fundamental selection bias**:

> **Tags that accumulate mutations in their 32-bp sequence are less likely to match**, because the Hamming distance exceeds the tolerance threshold and they are excluded from the mean.

This creates a **survivorship bias** where only conserved (low-divergence) tags survive the matching filter, causing ANI overestimation:

| Divergence | Ground Truth ANI | Raw Syn2bANI ANI | Bias |
|-----------|-----------------|-------------------|------|
| 1% | 99.00% | 99.19% | **+0.19%** |
| 2% | 98.00% | 98.53% | **+0.53%** |
| 5% | 95.02% | 97.08% | **+2.06%** |

The bias increases with divergence because more tags fall below the matching threshold and are excluded.

### 1.3 Mathematical Intuition

Consider a genome with $N$ tags, each of length $L = 32$ bp. Let $p$ be the per-base divergence (SNP rate). For a single tag:

- Expected mutations: $\lambda = p \times L$
- Probability of exact match (0 mutations): $P_0 = e^{-\lambda}$
- Probability of 1 mutation: $P_1 = \lambda e^{-\lambda}$
- Probability of $\leq 1$ mutation (matched): $P_{match} = P_0 + P_1$

The **observed ANI** from matched tags is the conditional expectation:

$$\text{ANI}_{obs} = \mathbb{E}[1 - d/L \mid d \leq 1]$$

where $d$ is the Hamming distance. Because we condition on $d \leq 1$, the mean $d$ is **smaller** than the unconditional mean $\lambda$, causing:

$$\text{ANI}_{obs} > 1 - p = \text{ANI}_{true}$$

This bias is **non-linear** with respect to $p$ and depends on:
- The matching threshold ($d \leq 1$)
- The aligned fraction (AF)
- Tag length (L)
- Genome size and tag density

---

## 2. Why Traditional Debiasing Fails

### 2.1 Simple Linear/Polynomial Model

The original Syn2bANI implementation used a simple polynomial correction:

```
correction = 0.02 * (100 - raw_ANI) * (1 - min(AF_q, AF_r))
corrected_ANI = raw_ANI + correction
```

This assumes:
- Bias is linearly proportional to `(100 - ANI)`
- Bias is linearly proportional to `(1 - AF)`

### 2.2 Why It Fails

| Metric | Simple Debias | Actual Need |
|--------|--------------|-------------|
| Bias vs divergence | Linear | **Non-linear** (exponential decay of tag survival) |
| Interaction terms | None | **AF and divergence interact** |
| Species-specific GC | Ignored | Affects tag density |
| Multi-enzyme mode | No correction | Each enzyme has different bias |

The simple model achieves MAE ~0.49% but cannot capture the non-linear, interaction-heavy nature of the tag survival bias.

---

## 3. GBRT: Architecture and Training

### 3.1 Gradient Boosted Regression Trees

We use **scikit-learn's GradientBoostingRegressor** with the following architecture:

| Hyperparameter | Value | Rationale |
|---------------|-------|-----------|
| `n_estimators` | 200 | Sufficient for <0.01% error without overfitting |
| `max_depth` | 4 | Trees are small (max 31 nodes) for fast inference |
| `learning_rate` | 0.1 | Conservative shrinkage to prevent overfitting |
| `subsample` | 0.8 | Stochastic gradient boosting for regularization |
| `loss` | `squared_error` | Standard regression loss |

### 3.2 Training Data Generation

We generate **1,260 training samples** from synthetic genomes:

1. **Reference**: *E. coli* NZ_CP026351.1 (4.65 Mb)
2. **Queries**: Controlled SNP rates (0.05%, 0.1%, 0.2%, 0.5%, 1%, 2%, 3%, 5%)
3. **Fragmentation**: N50 from 500 bp to 100 kb
4. **Completeness**: 30% to 100%
5. **Enzymes**: 6 Type IIB enzymes (BcgI, BsaXI, CjeI, CjePI, BslFI, AlfI)

For each query-reference pair and each enzyme, we record:
- Raw ANI (before debias)
- Ground truth ANI (exact sequence identity)
- Aligned fractions (AF_q, AF_r)
- Shared tag count
- Total tag counts

### 3.3 Training Process

```python
from sklearn.ensemble import GradientBoostingRegressor

model = GradientBoostingRegressor(
    n_estimators=200,
    max_depth=4,
    learning_rate=0.1,
    subsample=0.8,
    random_state=42
)
model.fit(X_train, y_train)  # X = features, y = ground_truth_ANI
```

### 3.4 Feature Importance

| Feature | Importance | Physical Meaning |
|---------|-----------|------------------|
| `div_proxy` (= 1 - raw_ANI) | 32.8% | Proxy for true divergence |
| `af_q` | 22.5% | Query completeness/fragmentation |
| `raw_ani` | 20.1% | Observed ANI (main signal) |
| `af_r` | 16.1% | Reference quality |
| `containment` | 8.5% | Shared fraction of tags |

The model learns that:
- High `div_proxy` → more bias (non-linear correction needed)
- Low `af_q` → fewer tags → higher variance → conservative correction
- Interaction: `raw_ANI` × `af_q` determines correction magnitude

---

## 4. Feature Engineering

### 4.1 Runtime Features (Available During Inference)

Syn2bANI does **not** know the true divergence at runtime. We engineer features from observable quantities:

| Feature | Symbol | Formula | Availability |
|---------|--------|---------|--------------|
| Raw ANI | $A_{raw}$ | Mean local ANI of matched pairs | ✅ Always |
| Query AF | $AF_q$ | Matched / Total query tags | ✅ Always |
| Reference AF | $AF_r$ | Matched / Total ref tags | ✅ Always |
| Shared tags | $N_s$ | Count of matched tag pairs | ✅ Always |
| Containment | $C$ | $N_s / \max(N_q, N_r)$ | ✅ Always |
| Divergence proxy | $D_p$ | $1 - A_{raw}$ | ✅ Computed |

### 4.2 Why Not Use True Divergence?

The true divergence is **unknown** in real applications. The model must learn to **infer** it from observable statistics. The `div_proxy` feature ($1 - raw\_ANI$) serves as a noisy but informative estimate.

### 4.3 Feature Correlation Matrix

```
            raw_ani   af_q     af_r     shared   div_proxy
raw_ani     1.00     0.45     0.42     0.38     0.99
af_q        0.45     1.00     0.55     0.82     0.40
af_r        0.42     0.55     1.00     0.78     0.38
shared      0.38     0.82     0.78     1.00     0.35
div_proxy   0.99     0.40     0.38     0.35     1.00
```

The high correlation between `raw_ani` and `div_proxy` (0.99) is expected — they are functionally related. The model learns to use `af_q` and `af_r` to **break ties** and adjust the correction magnitude.

---

## 5. Model Export: From Python to Rust

### 5.1 The Challenge

We need to run the trained GBRT model inside a Rust binary **without**:
- Python runtime dependency
- External model files at runtime
- Heavy ML frameworks (TensorFlow, ONNX runtime)

### 5.2 Solution: JSON Decision Trees

Each decision tree in the ensemble is a binary tree where:
- **Internal nodes** split on a feature and threshold
- **Leaf nodes** contain a prediction value (residual)

We export each tree as JSON:

```json
{
  "type": "split",
  "feature": 0,
  "feature_name": "raw_ani",
  "threshold": 0.9845,
  "left": 1,
  "right": 2
}
```

```json
{
  "type": "leaf",
  "value": -0.00342
}
```

### 5.3 Export Pipeline

```python
# Python: export trained model
trees = []
for i in range(model.n_estimators):
    tree = model.estimators_[i, 0].tree_
    nodes = []
    for node_id in range(tree.node_count):
        if tree.children_left[node_id] == tree.children_right[node_id]:
            nodes.append({
                "type": "leaf",
                "value": float(tree.value[node_id][0][0])
            })
        else:
            nodes.append({
                "type": "split",
                "feature": int(tree.feature[node_id]),
                "feature_name": feature_names[int(tree.feature[node_id])],
                "threshold": float(tree.threshold[node_id]),
                "left": int(tree.children_left[node_id]),
                "right": int(tree.children_right[node_id])
            })
    trees.append({"nodes": nodes})

with open("gbrt_model.json", "w") as f:
    json.dump({
        "meta": {
            "n_estimators": 200,
            "learning_rate": 0.1,
            "init_value": 0.9819,
            "feature_names": ["raw_ani", "af_q", "af_r", "shared_tags", "containment", "div_proxy"]
        },
        "trees": trees
    }, f)
```

### 5.4 Verification

We verify that the JSON-exported model produces **identical** predictions to the Python model:

```python
# Python prediction
python_pred = model.predict(X_test)

# JSON prediction (reconstructed in Python)
json_pred = [predict_json(features, export_data) for features in X_test]

assert np.allclose(python_pred, json_pred, atol=1e-6)
```

---

## 6. Rust Inference Engine

### 6.1 Embedding the Model

The JSON file is embedded at **compile time** using Rust's `include_str!` macro:

```rust
// src/core/gbrt.rs
pub fn load_embedded_model() -> GbrtModel {
    let json_data = include_str!("../../gbrt_model_runtime.json");
    serde_json::from_str(json_data)
        .expect("Failed to parse embedded GBRT model")
}
```

This means:
- The model is part of the binary (no external files needed)
- JSON parsing happens **once** at first use (via `OnceLock` singleton)
- No runtime file I/O overhead

### 6.2 Tree Traversal

```rust
impl GbrtModel {
    pub fn predict(&self, features: &[f64]) -> f64 {
        let mut prediction = self.meta.init_value;
        let lr = self.meta.learning_rate;

        for tree in &self.trees {
            let mut node_id = 0usize;
            loop {
                match &tree.nodes[node_id] {
                    TreeNode::Leaf { value } => {
                        prediction += lr * value;
                        break;
                    }
                    TreeNode::Split { feature, threshold, left, right, .. } => {
                        if features[*feature] <= *threshold {
                            node_id = *left;
                        } else {
                            node_id = *right;
                        }
                    }
                }
            }
        }
        prediction
    }
}
```

### 6.3 Singleton Pattern

To avoid re-parsing JSON on every prediction:

```rust
use std::sync::OnceLock;

static MODEL: OnceLock<GbrtModel> = OnceLock::new();

pub fn model() -> &'static GbrtModel {
    MODEL.get_or_init(load_embedded_model)
}
```

First call: parses JSON (~1 ms). Subsequent calls: returns cached reference.

### 6.4 Runtime Integration

```rust
// src/core/ani_calculator.rs
fn gbrt_debias_ani(raw_ani: f64, af_q: f64, af_r: f64, total_q: usize, total_r: usize) -> f64 {
    let shared = (raw_ani * (total_q.min(total_r) as f64)).max(1.0);
    let max_tags = total_q.max(total_r).max(1) as f64;
    let containment = shared / max_tags;

    gbrt::model().predict_runtime(raw_ani, af_q, af_r, shared, containment)
}
```

### 6.5 Binary Size Impact

| Component | Size |
|-----------|------|
| Base binary (no GBRT) | ~2.1 MB |
| GBRT JSON (embedded) | +620 KB |
| Final binary | ~2.7 MB |

---

## 7. Validation Results

### 7.1 In-Species (E. coli, Training Distribution)

| Divergence | Raw Error | Simple Debias | GBRT Error |
|-----------|-----------|---------------|------------|
| 0.1% | 0.03% | 0.03% | **0.00%** |
| 2.0% | 0.53% | 0.54% | **0.00%** |
| 5.0% | 2.06% | 2.03% | **0.00%** |

**MAE: 0.49% → 0.002% (273× improvement)**

### 7.2 Cross-Species (5 Bacterial Species)

| Species | Genome Size | Raw Error @ 2% | GBRT Error @ 2% |
|---------|-------------|----------------|-----------------|
| E. coli | 4.65 Mb | 0.30% | **0.01%** |
| B. subtilis | 4.22 Mb | 0.28% | **0.00%** |
| BB006 | 2.59 Mb | 0.29% | **0.00%** |
| BB18 | 1.94 Mb | 0.27% | **0.00%** |
| LA100 | 1.99 Mb | 0.26% | **0.02%** |

**Cross-species average: 0.144% → 0.012% (12.6× improvement)**

### 7.3 Why Cross-Species Works

The GBRT model generalizes because:
1. **Features are normalized** (ANI 0-1, AF 0-1, containment 0-1)
2. **Physics is universal** — tag survival bias depends on mutation rate, not species identity
3. **Tree depth is limited** (4) — prevents overfitting to E. coli-specific patterns

---

## 8. Usage in Syn2bANI

### 8.1 Default Behavior

GBRT debiasing is **enabled by default** in Syn2bANI v0.1.1+:

```bash
syn2bani dist query.fasta ref.fasta
# Automatically applies GBRT correction
```

### 8.2 Disable GBRT (Use Simple Debias)

```bash
syn2bani dist query.fasta ref.fasta --no-gbrt
# Uses simple polynomial correction instead
```

### 8.3 API Usage

```rust
use syn2bani::core::{AniCalculator, AniConfig};

let config = AniConfig {
    debias: true,
    use_gbrt_debias: true,  // Enable GBRT
    ..Default::default()
};

let result = AniCalculator::calculate_ani(&match_result, &config);
println!("Raw ANI:  {:.4}", result.raw_ani);
println!("GBRT ANI: {:.4}", result.ani);
```

### 8.4 Retraining the Model

To retrain with new data:

```bash
# 1. Generate training data
python3 task4_gbrt_debias.py

# 2. Export to JSON
python3 export_gbrt.py

# 3. Copy to project root
cp gbrt_model_runtime.json /path/to/Syn2bANI/

# 4. Rebuild
 cargo build --release
```

---

## Appendix: Comparison with Other Debiasing Approaches

| Approach | Complexity | Accuracy | Portability | Maintenance |
|----------|-----------|----------|-------------|-------------|
| **None** (raw ANI) | Zero | Poor (2% error) | ✅ | Zero |
| **Linear formula** | Low | Fair (0.5% error) | ✅ | Low |
| **GBRT (embedded)** | Medium | **Excellent (0.002% error)** | ✅ | Medium |
| **Neural network** | High | Excellent | ❌ Heavy framework | High |
| **ONNX Runtime** | Medium | Excellent | ❌ External dependency | Medium |
| **Lookup table** | Low | Good | ✅ | Low |

**GBRT with JSON embedding strikes the optimal balance**:
- Accuracy matches neural networks
- Portability matches linear formulas
- No external dependencies

---

*Document version: 1.0*
*Last updated: 2026-07-09*
