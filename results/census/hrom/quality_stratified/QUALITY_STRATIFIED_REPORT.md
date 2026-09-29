# HROM CheckM/quality sensitivity analysis

- pairs seen: 62,872,901
- qualifying pairs: 60,014,738
- top 5% pairs per axis: 3,000,736

## Quality-stratified mode counts

| stratum | pairs | species tested | structural | ANI | both | not |
|---|---:|---:|---:|---:|---:|---:|
| all_qualified | 60,014,738 | 521 | 96 | 161 | 36 | 228 |
| high_quality_comp90_cont5 | 18,210,578 | 354 | 60 | 93 | 19 | 182 |
| medium_comp70_cont10 | 43,676,804 | 481 | 88 | 139 | 27 | 227 |
| high_completeness_only | 18,210,578 | 354 | 60 | 93 | 19 | 182 |
| continuous_N50_ge10kb | 23,015,221 | 425 | 64 | 118 | 20 | 223 |

## High-quality mode stability

| mode | all | high-quality | intersection | Jaccard |
|---|---:|---:|---:|---:|
| structural_enriched | 96 | 60 | 53 | 0.515 |
| ani_enriched | 161 | 93 | 79 | 0.451 |
| both_enriched | 36 | 19 | 9 | 0.196 |

See TSV files for exact odds ratios, P values, per-stratum species and robust candidates.
