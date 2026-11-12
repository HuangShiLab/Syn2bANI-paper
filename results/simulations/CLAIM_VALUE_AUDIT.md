# Exact-truth simulation claim-value audit, 2026-11-12

This audit reconciles the manuscript and current Figure 2 with the packaged
source tables. The current comparator tier is skani 0.3.2 with default
parameters, matching Online Methods.

## Twelve-pair ANI ladder

Truth is encoded in `simindel/manifest.tsv`; Syn2bANI estimates are in
`simindel_results_4e.tsv`; comparator calls are in
`simindel_cross_tool_4e.tsv`.

| method | n | MAE (ANI points) | mean signed error |
|---|---:|---:|---:|
| Syn2bANI | 12 | 0.0732 | +0.0647 |
| skani 0.3.2 | 12 | 0.4775 | -0.4775 |
| FastANI | 12 | 0.7423 | -0.7423 |

The earlier skani 0.1.0 result (MAE 0.3775) is retained in
`SKANI_0.3.2_RERUN.md` as a historical comparison and is no longer quoted as
the manuscript headline.

## Robustness families

| family | source | n | MAE | signed-error range |
|---|---|---:|---:|---|
| indel sweep | `simindel_sweep_results_4e.tsv` | 5 | 0.0807 | +0.0322 to +0.1844 |
| fragmentation | `simfrag_fixed.tsv` | 8 rows; four distinct contig levels | 0.0926 | +0.0176 to +0.1794 |
| accessory content | `simacc_fixed.tsv` | 6 | 0.1121 | +0.0382 to +0.2588 |

The manuscript rounds these to 0.081, 0.093 and 0.114. Fragmentation rows are
not deduplicated before plotting because forward and shuffled controls have
identical estimates; grouping by contig count gives the same four error values
and the same 0.093 MAE as the manuscript. The older
`results/SIMULATION_RERUN_REPORT.md` value 0.0695 used the duplicated row set
directly and is superseded for manuscript reporting.

## Figure and caption

The packaged figure is
`paper/figures/main/final/figure2.pdf`; the current rendered labels report
Syn2bANI 0.073, skani 0.477 and FastANI 0.742. The manuscript Abstract,
Results and Figure 2 caption were updated to the same values on 2026-11-12.

## Additional release audit

Mosaic simulation outputs and their manifest are now packaged in
`results/simulations/` for Supplementary Figure S6. The S. aureus supplement
figure was regenerated from `case_studies/fda_argos_s_aureus/results/triangle_rerun.tsv`
and reports the current Syn2bANI v0.1.1 result: 2,211 pairs, 1,337/2,211
(60.47%) at zero breakpoints, maximum seven breakpoints, and ANI 93.6-100.0%.
The obsolete generator using the pre-v0.1.1 triangle was replaced.
