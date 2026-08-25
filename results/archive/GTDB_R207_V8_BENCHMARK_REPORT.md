# GTDB-R207 v8 Benchmark Report

Matrix: `results/matrix_gtdb_r207_100k_v8_final.tsv`
Total pairs: 45000
v8 reported pairs: 9739
skani reported pairs: 4467
Pairs with FastANI reference: 672

## Overall metrics vs FastANI

| Method | n | MAE | RMSE | Pearson r | bias |
|--------|---|-----|------|-----------|------|
| Syn2bANI v8 gamma | 620 | 3.304% | 4.249% | 0.7262 | +3.018% |
| Syn2bANI v8 uniform | 620 | 4.235% | 4.931% | 0.8079 | +4.235% |
| skani | 672 | 0.608% | 0.914% | 0.9798 | +0.295% |

## Metrics by flag (v8 gamma vs FastANI)

| flag | n | MAE | RMSE | r | bias |
|------|---|-----|------|---|------|
| ok | 215 | 4.862% | 5.690% | 0.7787 | +4.862% |
| INCONSISTENT | 405 | 2.476% | 3.232% | 0.7907 | +2.039% |

## Flag counts

| flag | count |
|------|-------|
| BELOW_DETECTION | 35261 |
| ok | 5252 |
| INCONSISTENT | 4487 |

## Per-phylum MAE (>=10 pairs with FastANI ref)

| phylum | n | MAE | bias |
|--------|---|-----|------|
| p__Actinobacteriota | 62 | 2.720% | +2.375% |
| p__Bacteroidota | 36 | 2.788% | +2.338% |
| p__Chloroflexota | 12 | 4.516% | +4.516% |
| p__Cyanobacteria | 27 | 3.744% | +3.706% |
| p__Desulfobacterota | 20 | 4.313% | +4.313% |
| p__Firmicutes | 39 | 2.888% | +2.317% |
| p__Firmicutes_A | 53 | 3.524% | +3.421% |
| p__Halobacteriota | 18 | 2.577% | +2.322% |
| p__Myxococcota | 10 | 3.839% | +3.839% |
| p__Planctomycetota | 10 | 3.666% | +3.063% |
| p__Proteobacteria | 191 | 2.864% | +2.422% |
| p__Thermoproteota | 10 | 4.013% | +4.013% |
| p__Verrucomicrobiota | 11 | 4.643% | +4.643% |

## Interpretation

- v8 gamma shows a systematic +3.0% bias vs FastANI across GTDB-R207 mid-high pairs.
- The bias is present in all major phyla (2.4–4.6%), suggesting a global offset rather than GC-specific effect.
- skani agrees much better with FastANI (MAE 0.61%), but both are k-mer methods and may share similar biases at this ANI range.
- Notably, "ok" pairs have larger MAE (4.86%) than "INCONSISTENT" pairs (2.48%), opposite to the oral/gut validation.
- This discrepancy vs simulations (0.06–0.36% MAE on exact ground truth) strongly suggests the issue is real-draft assembly effects or the FastANI reference itself, not the core MLE model.
- Independent ground truth (ANIm/nucmer or minimap2) is required to determine whether v8 or FastANI is closer to truth.