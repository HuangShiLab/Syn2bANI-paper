# SV re-analysis — Syn2bANI v0.1.1 breakpoint_count (43,334 pairs)

- breakpoint_count: median 2, mean 7.1, max 269
- correlation with contig count: r = -0.049, rho = -0.028 (v0.1.1 counter)
- synteny_blocks vs contig count: r = 0.534 (new primary-block definition)

| target | raw r | \|ANIm | \|n_contigs | \|both | n |
|---|---:|---:|---:|---:|---:|
| dnadiff breakpoints (all) | 0.361 | 0.369 | 0.398 | **0.411** | 43,334 |
| dnadiff breakpoints (>=10 kb) | 0.221 | 0.229 | 0.263 | **0.276** | 43,334 |
| dnadiff large indels (>=10 kb) | 0.289 | 0.300 | 0.313 | **0.327** | 43,334 |
| dnadiff inversions (all) | 0.495 | 0.499 | 0.495 | **0.500** | 43,334 |

- >=95% ANIm subset (n = 654): Spearman rho vs large indels = 0.633 (withdrawn value was 0.674)
- cross-implementation concordance on the same 43,312 pairs: Spearman rho = 0.863, Pearson r = 0.843; median syn2bani/syn2b count ratio = 0.30

Replaces the withdrawn r = 0.414/0.453 table in SV_REANALYSIS.md.
