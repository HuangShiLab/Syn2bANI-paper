# SV method comparison on GTDB 50k held-out pairs
Pairs with Syn2bANI output: 43,334

## Breakpoint count correlations with Syn2bANI breakpoint_count
- dnadiff (all gaps): Pearson r=0.4649, Spearman r=0.5271, MAE=424.9
- dnadiff synteny (all gaps): Pearson r=-0.1360, Spearman r=-0.1238, MAE=30.5
- minimap2: Pearson r=0.3004, Spearman r=0.3566, MAE=179.5
- minimap2 synteny: Pearson r=-0.0191, Spearman r=0.1902, MAE=30.6
- dnadiff min-gap 10000 bp: Pearson r=0.4894, Spearman r=0.5550, MAE=401.4
- dnadiff large indels min-gap 10000 bp: Pearson r=0.5522, Spearman r=0.6946, MAE=255.9
- dnadiff min-gap 5000 bp: Pearson r=0.4793, Spearman r=0.5438, MAE=412.1
- dnadiff large indels min-gap 5000 bp: Pearson r=0.5087, Spearman r=0.6588, MAE=326.6

## Correlations with alignment-based synteny/coverage scores

### Syn2bANI synteny_score (anchor-adjacency conservation)
- dnadiff synteny (all gaps): Pearson r=0.1193, Spearman r=0.3292
- minimap2 synteny: Pearson r=0.0181, Spearman r=-0.0129

### Syn2bANI af_query (base-pair chain coverage)
- dnadiff synteny (all gaps): Pearson r=0.3539, Spearman r=0.3939
- minimap2 synteny: Pearson r=0.2305, Spearman r=0.2627

### Syn2bANI synteny_blocks
- dnadiff synteny (all gaps): Pearson r=-0.1680, Spearman r=-0.2133
- minimap2 synteny: Pearson r=-0.1674, Spearman r=-0.0530

## Summary statistics
- breakpoint_count: mean=30.9, median=20.0
- dnadiff_breakpoints: mean=455.7, median=333.0
- mm2_breakpoints: mean=209.3, median=156.0
