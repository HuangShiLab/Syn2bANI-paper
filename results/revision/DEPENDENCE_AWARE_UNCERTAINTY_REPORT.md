# Dependence-aware uncertainty (revision audit)

## Scope and interpretation

GTDB benchmark intervals use 2,000 genome-level bootstrap replicates: whole genomes are resampled with replacement and pairs inherit the product of endpoint multiplicities, so all observations containing a genome move together. HROM and OAPGC intervals use 2,000 species-level bootstrap replicates. Species-level resampling quantifies composition uncertainty but does not remove pair-level dependence, shared ancestry, assembly-quality correlation or database sampling structure; these remain explicit limitations.

## GTDB-R207 43,334-pair benchmark

MAE differences are reported on the common complete-case subset (39,903 pairs). Positive values mean the second method has lower MAE. The point estimates are reproducible from `results/gtdb50k/s2b_50k.tsv` with v6 calibration (`ani_cal`) and the operational hybrid rule (raw gated at ANI estimate ≥98%).

- `mae_syn2bani_hybrid_minus_skani` = -0.337 (genome-bootstrap 95% CI -0.359 to -0.315)
- `mae_syn2bani_v6_minus_skani` = -0.338 (genome-bootstrap 95% CI -0.360 to -0.317)

## OAPGC oral-airway site test

The sequence-distance metric is `100 - skani ANI`; it is an ANI-derived distance rather than a base-level SNP count.

- `median_sequence_distance_auc` = 0.574 (species-bootstrap 95% CI 0.555 to 0.582)
- `median_BP_auc` = 0.494 (species-bootstrap 95% CI 0.489 to 0.509)

## HROM candidate divergence-mode fractions

- `fraction_ani_enriched` = 0.307 (species-bootstrap 95% CI 0.269 to 0.346)
- `fraction_both_enriched` = 0.069 (species-bootstrap 95% CI 0.048 to 0.092)
- `fraction_not_enriched` = 0.440 (species-bootstrap 95% CI 0.395 to 0.486)
- `fraction_structural_enriched` = 0.184 (species-bootstrap 95% CI 0.154 to 0.217)

Machine-readable intervals: `results/revision/dependence_aware_uncertainty.tsv`.
