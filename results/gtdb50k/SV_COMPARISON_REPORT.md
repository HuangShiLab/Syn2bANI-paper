# SV method comparison on GTDB 50k held-out pairs
Pairs with Syn2bANI output: 43,334

## Breakpoint count correlations with Syn2bANI breakpoint_count
- dnadiff (all gaps): Pearson r=0.1327, Spearman r=0.2948, MAE=434.6
- minimap2: Pearson r=0.0641, Spearman r=0.1762, MAE=190.9
- minimap2 synteny: Pearson r=0.1484, Spearman r=0.4350, MAE=23.4
- dnadiff min-gap 10000 bp: Pearson r=0.1469, Spearman r=0.3187, MAE=411.2
- dnadiff large indels min-gap 10000 bp: Pearson r=0.2553, Spearman r=0.5628, MAE=264.6
- dnadiff min-gap 5000 bp: Pearson r=0.1408, Spearman r=0.3086, MAE=421.8
- dnadiff large indels min-gap 5000 bp: Pearson r=0.3178, Spearman r=0.6276, MAE=334.7

## Correlations with alignment-based synteny/coverage scores

### Syn2bANI anchor_adjacency (anchor-adjacency conservation)

### Syn2bANI af_query (base-pair chain coverage)

### Syn2bANI synteny_blocks

## Summary statistics
- breakpoint_count: mean=23.6, median=12.0
- dnadiff_breakpoints: mean=455.7, median=333.0
- mm2_breakpoints: mean=209.3, median=156.0
