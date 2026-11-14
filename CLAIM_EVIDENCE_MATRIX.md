# Syn2bANI claim-evidence matrix

Status legend: **pass** means current text, source table and display agree after
the 2026-11-12 skill audit.

| ID | Claim | Scope/denominator | Expected direction | Source | Display | Status |
|---|---|---|---|---|---|---|
| C01 | Exact-truth MAE 0.073 | 12 deterministic ladder pairs | Syn2bANI lower MAE | `results/simulations/simindel_results_4e.tsv` + manifest | Fig. 2a,b | pass |
| C02 | skani 0.3.2 exact-truth MAE 0.477 | same 12 pairs | comparator higher MAE | `results/simulations/simindel_cross_tool_4e.tsv` | Fig. 2a,b | fixed 2026-11-12 |
| C03 | FastANI exact-truth MAE 0.742 | same 12 pairs | comparator higher MAE | same as C02 | Fig. 2a,b | pass |
| C04 | Robustness MAE 0.081/0.093/0.114 | indel/fragmentation/accessory sweeps | low error; AF tracks shared fraction | `results/simulations/*fixed.tsv` | Fig. 2c-e | pass |
| C05 | GTDB held-out MAE 0.619 versus 0.957/0.977 | 39,903 complete-case pairs of 43,334 | calibrated Syn2bANI lower MAE | `results/gtdb50k/GTDB50K_HELDOUT_REPORT.md` | Fig. 3; Table 1 | pass |
| C06 | Comparator calibration narrows/reverses advantage | common complete-case pairs | gain mostly bias removal | `results/gtdb50k/CALIBRATED_COMPARATORS.md` | Table 1 | pass |
| C07 | Hybrid 0.615 on unified benchmark | 40,629 complete-case pairs | high-ANI recovery without low-ANI loss | `results/gtdb50k/UNIFIED_BENCHMARK_REPORT.md` | Fig. 3; Table 1 | pass |
| C08 | Genome-bootstrap v6-minus-skani −0.338 | 39,903 pairs, genome-level resampling | negative interval | `results/revision/dependence_aware_uncertainty.tsv` | Results | pass |
| C09 | SV follow-up ~60–90× faster | 7 top-hit pairs | Syn2bANI faster | `results/one_to_all/sv_benchmark.tsv` | ED Fig. 1f; Table 2 | pass |
| C10 | Breakpoint partial correlation 0.411 | 43,334 held-out pairs, controlling ANIm and contigs | moderate not exact | `results/gtdb50k/rerun_v011/s2b_50k_v011_merged.tsv` | Fig. 4c | pass |
| C11 | High-ANI inversion agreement r=0.996 | 3,826 pairs at ANIm ≥97% | near identity | `results/gtdb50k/inverted_fraction_comparison_report.md` | Fig. 4b | regression direction clarified |
| C12 | abfA periphery coverage is lower for curated deleted/peripheral-truncated strains | 183 non-reference genomes: 127 complete and 56 deleted/truncated | deleted lower coverage | `results/b_longum_abfA/abfA_region_coverage.tsv`; `results/b_longum_abfA/ABFA_CASE_REPORT.md` | Fig. 5 | pass |
| C13 | FDA-ARGOS 98% structural and 58% genotype discordance | 103 pairs at ANI ≥99.9 | high ANI does not imply equivalence | `case_studies/fda_argos/validate/VALIDATION.md`, `amr_discordance.tsv` | Results; Supp Fig. S15 | pass |
| C14 | cagPAI lineage association survives; disease association does not | 528 isolates; lineage-stratified tests | lineage confounding explicit | `case_studies/h_pylori_cagpai/results/` | Fig. 6 | pass |
| C15 | HROM enrichment fractions | 60,014,738 qualifying pairs; 521 tested clusters | axes provide partly distinct information | `results/census/hrom/hypermodes/` | ED Fig. 3 | pass |
| C16 | OAPGC sequence-distance enrichment without breakpoint excess | 200 species; 1,264,327 cross-site and 16,067,030 within-habitat pairs | ANI-distance positive; breakpoints null | `results/census/oapgc/site_effects/` | ED Fig. 4 | pass |
| C17 | GTDB catalog-scale decline at high ANI | 1,311 complete species with ≥100 expected pairs | breakpoint-positive fraction declines at ≥99% | `results/census/gtdb/gtdb_species_ani_sv_joint.tsv` | ED Fig. 5 | pass |
