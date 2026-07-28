# Oral/Gut Independent Validation Report (Syn2bANI v8) — flag-based analysis

Total pairs: 1225
Reported pairs (ok + INCONSISTENT): 122
Pairs with FastANI reference: 100

## Metrics vs FastANI by flag

| flag | n | model | MAE | RMSE | r | bias |
|------|---|-------|-----|------|---|------|
| all reported | 100 | gamma | 0.552% | 1.056% | 0.9804 | +0.421% |
| all reported | 100 | uniform | 1.165% | 2.019% | 0.9745 | +1.165% |
| ok | 66 | gamma | 0.293% | 0.343% | 0.9034 | +0.240% |
| ok | 66 | uniform | 0.553% | 0.617% | 0.9396 | +0.553% |
| INCONSISTENT | 34 | gamma | 1.054% | 1.747% | 0.9699 | +0.773% |
| INCONSISTENT | 34 | uniform | 2.351% | 3.354% | 0.9585 | +2.351% |

## Key findings

- **Diagnostic flag works**: ok-only MAE 0.293% vs all-reported 0.552%.
- INCONSISTENT pairs have ~3.6× higher MAE than ok pairs, confirming the built-in QC identifies unreliable estimates on real data.
- retention vs |err|: r=-0.6880 — lower retention predicts larger error.
- No high pair had het_shape = uniform (LRT always accepted heterogeneity), so the residual +0.240% bias in ok pairs is not from overly conservative gating.
- The +0.421% overall bias is driven by INCONSISTENT outliers; ok-only bias is +0.240%.

## Caveats
- Only 100 high-ANI (same-species) pairs have FastANI reference.
- 1,100 low-ANI pairs lack a reference; FastANI does not report below ~80% ANI.
- True mid/low-ANI validation requires ANIm/nucmer or minimap2 ground truth.