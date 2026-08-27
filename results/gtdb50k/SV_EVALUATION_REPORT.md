# SV truth evaluation against dnadiff 1-to-1 coords (GTDB 50k held-out)
Pairs with parseable dd.1coords: 43,334
Pairs with syn2bani output: 43,334

## Overall metrics
- breakpoint_count MAE vs dnadiff: 424.927
- breakpoint_count Pearson r: 0.465
- anchor_adjacency Pearson r: 0.119
- Rearrangement detection (truth > 0 vs pred > 0): precision=0.999, recall=0.931, F1=0.964, specificity=0.608
  TP=40269, FP=40, FN=2963, TN=62

## Per-band breakpoint_count MAE
- 80-85: 461.413 (n=12172)
- 85-90: 496.363 (n=16000)
- 90-95: 322.427 (n=14758)
- 95-100: 240.782 (n=404)
