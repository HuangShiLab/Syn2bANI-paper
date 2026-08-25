# GBRT v5 Enzyme Comparison on Mid-ANI Independent Validation

## Setup
- Models trained on GTDB-R207 100k pair matrix (`train_bcgi.tsv`, `train_cjepi.tsv`, `train_combined.tsv`)
- Training samples: 728 pairs with FastANI reference
- Validation: oral/gut independent pairs (`val_combined.tsv`) filtered to 15 pairs where all tools returned results
- SLURM evaluation job: `3807893` (completed)
- Evaluation script: `scripts/evaluate_gbrt_v5_combined.py`

## Training in-sample performance

| Model | Best params | In-sample MAE vs FastANI | In-sample RMSE |
|-------|-------------|--------------------------|----------------|
| BcgI  | 1000 trees, depth 6, lr 0.02 | 0.278% | 0.763% |
| CjePI | 1000 trees, depth 5, lr 0.02 | 0.261% | 0.324% |
| Combined | 300 trees, depth 5, lr 0.05 | 0.227% | 0.284% |

## Independent validation performance (n=15 pairs)

| Method | MAE vs FastANI | RMSE vs FastANI |
|--------|----------------|-----------------|
| Syn2bANI BcgI raw | 10.051% | 10.054% |
| Syn2bANI BcgI mash | 6.656% | 6.753% |
| **Syn2bANI GBRT v5 BcgI** | **3.107%** | **3.403%** |
| Syn2bANI CjePI raw | 9.587% | 9.593% |
| Syn2bANI CjePI mash | 6.507% | 6.556% |
| Syn2bANI GBRT v5 CjePI | 3.699% | 4.096% |
| Syn2bANI GBRT v5 Combined | 3.417% | 3.764% |
| **skani** | **0.468%** | **0.524%** |

## Key observations
1. **BcgI gives the best corrected result** (3.11% MAE), slightly better than combined (3.42%) and CjePI (3.70%).
2. **Mash-like ANI already reduces the raw overestimation** by ~3.4 percentage points (BcgI) and ~3.1 percentage points (CjePI).
3. **GBRT correction further halves the error**, but the residual error is still ~6-7x larger than skani.
4. The validation set is **extremely small (15 pairs)** and **taxonomically narrow**: only *Bifidobacterium* and *Veillonella* pairs, all in the 86.5-87.6% FastANI range.
5. **Domain shift**: the training set contains zero *Bifidobacterium* or *Veillonella* pairs. The poor generalization is therefore partly a taxonomic domain-shift problem, not purely algorithmic.
6. **CjePI `mash_ani` is highly informative in training** (38.2% importance) but the CjePI-only model still generalizes worse than BcgI on this validation set, suggesting the extra information does not transfer well to these unseen genera.
7. The **combined model did not meaningfully use CjePI features** (combined feature importances: BcgI raw 44.9%, BcgI shared_log 52.7%, all CjePI features combined <3%).

## Interpretation
The large gap between in-sample MAE (~0.25%) and independent validation MAE (~3.1-3.7%) indicates that the current GBRT v5 models are **overfitting to the taxa and feature distribution of the GTDB-R207 training slice**, rather than learning a robust, genus-agnostic correction.

## Recommendations
1. **Do not embed the combined GBRT v5 model into Rust yet.** The validation error is too high and the model does not clearly outperform the single-enzyme BcgI model.
2. **Expand the independent validation set** to include more genera and more pairs in the 85-95% ANI range, ideally including taxa represented in training.
3. **Re-examine the training strategy**:
   - Ensure training covers the low-ANI tail (85-88%) with sufficient pairs from diverse genera.
   - Add regularization (shallower trees, more subsampling) or try a simpler linear/elastic-net model to reduce overfitting.
4. **Consider whether the fundamental issue is tag density/enzyme choice** before adding more training data: BcgI performs best, but the raw overestimation is still large. The Mash-like estimator helps, but a genus-aware correction may be needed.
5. If the final goal is to beat or match skani, the correction must be **much more conservative** and validated on a larger, taxonomically diverse hold-out set.
