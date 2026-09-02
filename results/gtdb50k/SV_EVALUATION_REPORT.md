# SV truth evaluation against dnadiff 1-to-1 coords (GTDB 50k held-out)
Pairs with parseable dd.1coords: 43,334
Pairs with syn2bani output: 43,334

## Overall metrics
- breakpoint_count MAE vs dnadiff: 434.558
- breakpoint_count Pearson r: 0.133
- anchor_adjacency Pearson r: 0.119
- Rearrangement detection (truth > 0 vs pred > 0): precision=0.999, recall=0.809, F1=0.894, specificity=0.608
  TP=34982, FP=40, FN=8250, TN=62

## Per-band breakpoint_count MAE
- 80-85: 466.528 (n=12172)
- 85-90: 505.021 (n=16000)
- 90-95: 336.444 (n=14758)
- 95-100: 264.824 (n=404)
