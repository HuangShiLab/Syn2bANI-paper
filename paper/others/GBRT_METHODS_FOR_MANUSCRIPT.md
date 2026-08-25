# Methods: GBRT-Based ANI Debiasing in Syn2bANI

> Condensed version suitable for the Methods section of a manuscript.

## Overview of the Fixed-Anchor Tag Matching Bias

Syn2bANI estimates ANI by comparing 2bRAD tags (~32 bp DNA fragments flanking Type IIB restriction sites) between query and reference genomes. The raw ANI is computed as the mean local identity of all matched tag pairs, where each local identity is defined as $1 - d/32$ with $d$ being the Hamming distance between the two tags.

This approach introduces a **systematic positive bias** because tags that accumulate mutations ($d \geq 2$) fail to match and are excluded from the mean. This creates a survivorship bias where only conserved tags contribute to the ANI estimate, causing overestimation that increases with true divergence (Table 1).

**Table 1. Systematic bias in raw Syn2bANI ANI estimates.**

| True divergence | Ground truth ANI | Raw Syn2bANI ANI | Bias |
|-----------------|------------------|-------------------|------|
| 0.1% | 99.90% | 99.93% | +0.03% |
| 1.0% | 99.00% | 99.19% | +0.19% |
| 2.0% | 98.00% | 98.53% | +0.53% |
| 5.0% | 95.02% | 97.08% | +2.06% |

For a tag of length $L = 32$ and per-base divergence $p$, the probability of a tag surviving the matching threshold ($d \leq 1$) is:

$$P_{match} = e^{-32p}(1 + 32p)$$

The observed ANI is the conditional expectation of identity given survival, which is always greater than the true identity $1-p$ because the conditioning excludes high-divergence tags.

## GBRT Debiasing Model

To correct this systematic bias, we trained a **Gradient Boosted Regression Tree (GBRT)** model to predict the true ANI from observable statistics.

### Model Architecture

We used scikit-learn's `GradientBoostingRegressor` with the following hyperparameters:
- **n_estimators**: 200 (number of trees)
- **max_depth**: 4 (maximum depth per tree, 31 nodes maximum)
- **learning_rate**: 0.1 (shrinkage factor)
- **subsample**: 0.8 (stochastic gradient boosting for regularization)
- **loss function**: squared error

### Training Data

Training data was generated from synthetic bacterial genomes derived from *E. coli* K-12 (NZ_CP026351.1, 4.65 Mb). We introduced controlled SNPs at rates ranging from 0.05% to 5%, fragmented the genomes at N50 values from 500 bp to 100 kb, and varied completeness from 30% to 100%. For each condition, we ran Syn2bANI with 6 different Type IIB enzymes (BcgI, BsaXI, CjeI, CjePI, BslFI, AlfI), yielding **1,260 training samples**. Each sample was labeled with the exact ground truth ANI (position-by-position sequence identity).

### Features

The model uses 6 features that are all available at runtime without knowing the true divergence:

| Feature | Description |
|---------|-------------|
| `raw_ani` | Observed ANI from matched tag pairs (0–1) |
| `af_q` | Aligned fraction of query (matched tags / total query tags) |
| `af_r` | Aligned fraction of reference (matched tags / total reference tags) |
| `shared_tags` | Number of matched tag pairs |
| `containment` | Shared tags / max(total query tags, total reference tags) |
| `div_proxy` | $1 - raw\_ani$ (proxy for true divergence) |

### Model Export and Embedding

The trained Python model was exported as **JSON-encoded decision trees** and embedded directly into the Rust binary using the `include_str!` compile-time macro. Each tree was serialized as a sequence of nodes with `split` (feature, threshold, left child, right child) or `leaf` (value) types. This approach eliminates all runtime dependencies on Python or machine learning frameworks.

### Inference in Rust

At runtime, the embedded JSON is parsed once into a `GbrtModel` singleton (via `std::sync::OnceLock`). Prediction is performed by traversing each of the 200 trees: starting at the root node, following the `feature <= threshold` branch until reaching a leaf, and accumulating the leaf value scaled by the learning rate. The final prediction is the sum of the initial value and all tree contributions. The total inference time is <1 μs per prediction.

### Model Validation

The model was validated on held-out test data from the same species (in-species validation) and on 4 additional bacterial species (cross-species validation). The GBRT model reduced the mean absolute error (MAE) from 0.49% (simple polynomial debias) to **0.002%** on the training distribution and **0.012%** across 5 species, demonstrating excellent generalization.

---

## Figure Legend (for manuscript)

**Figure: GBRT Debiasing Corrects Systematic ANI Overestimation.**
*(a)* Comparison of raw ANI, simple polynomial debias, and GBRT-corrected ANI against ground truth across 0.05%–5% divergence. *(b)* Absolute error for each method. *(c)* Cross-species validation on 5 bacterial genomes (E. coli, B. subtilis, and 3 environmental isolates). *(d)* Feature importance of the GBRT model showing that divergence proxy (1 - raw_ANI) and query AF are the most informative features.
