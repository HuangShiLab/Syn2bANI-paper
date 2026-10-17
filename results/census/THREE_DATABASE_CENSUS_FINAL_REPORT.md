# Three-database within-species census: final report

## Scope

This report integrates the final species-level structural/ANI censuses for GTDB-R207, HROM and OAPGC. All pair counts are unique unordered genome pairs. GTDB sensitivity uses both the all-reported denominator and the subset with numeric breakpoint values.

## GTDB-R207

- Effective scope: 261,737 genomes and 22,535 species clusters.
- Planned effective pairs: 653,584,719.
- Reported structural pairs after canonicalization/deduplication: 651,961,758 (99.7517%).
- Pairs with at least two breakpoints: 508,066,874 (77.9289% of reported pairs).
- Pairs without a numeric breakpoint value: 65,698,130 (10.0771%); breakpoint-positive fraction among breakpoint-informative pairs: 86.6618%.
- ANI pairs: 653,522,280 (99.9904% of the effective scope), all within expected species; no duplicate, self or cross-species pairs; no invalid ANI values.
- Final joint scope: 21,973 usable species clusters — 21,968 complete, 4 ANI-complete/structural-partial and 1 partial on both axes.
- The only structural shortfall is *Escherichia coli*: 337,170,504/338,793,465 pairs (99.5210%). Final block diagnostics found no fully missing block pair, no extra rows and no invalid identifiers; 1,407/1,431 block pairs were complete and 24 were partial. The residual deficit is therefore reported explicitly rather than imputed.
- Among 1,311 complete species with at least 100 expected pairs, median structural-positive fraction was 61.34% and median mean skani ANI was 98.47%. The median structural-positive fraction was 61.33% for 1,309 species at mean ANI at least 95%, 38.53% for 391 species at ANI at least 99%, and 9.52% for 53 species at ANI at least 99.9%.
- Quality sensitivity: 645,271,663 pairs remained at completeness at least 90% and contamination at most 5%; the breakpoint-positive fraction was 77.9002% versus 77.9289% for all reported pairs. Minimum pair N50 at least 10 kb gave 78.0244%, and complete HQ/N50 at least 10 kb gave 77.9773%.

## HROM

- Analysed genomes: 142,277; multi-genome species clusters: 2,241.
- Unique pairs with both structural and ANI measurements: 62,872,901.
- Pairs with at least two breakpoints: 17,535,013 (27.89%).
- Among 521 species with at least 200 qualifying pairs: 96 structural-similarity-enriched, 160 ANI-similarity-enriched, 36 both-similarity-enriched and 229 not enriched.
- High-quality subset (completeness at least 90%, contamination at most 5%): 18,210,578 pairs and 354 species; mode counts were 60 structural, 93 ANI, 19 both and 182 not enriched. Moderate overlap supports reporting both all-pair and high-quality lists rather than selecting the more biologically attractive one.

## OAPGC

- Analysed genomes: 98,229; multi-genome SGB clusters: 1,488.
- Informative structural pairs: 30,559,312; pairs with at least two breakpoints: 15,247,934 (49.90%).
- Oral-airway test: 1,264,327 oral-airway pairs versus 16,067,030 within-site pairs across 200 eligible species.
- Median species AUC: SNP distance 0.574; breakpoint count 0.494.
- Direction counts: SNP distance greater than 0.5 in 143/200 species (sign-test p = 1.0e-9); breakpoint count greater than 0.5 in 93/200 (p = 0.358).
- At BH q < 0.05: 73 SNP-distance-enriched only, 19 breakpoint-enriched only, 46 both and 62 not significant/depleted.
- The recent-strain subset (ANI at least 99.9%) had only 49 eligible oral-airway pairs and is a power audit, not evidence for recent transmission or directionality.

## Paper allocation

1. **Syn2bANI (Nature Methods target).** Use the three-database census as method-scale and discovery evidence: GTDB for hundreds of millions of pairs and quality sensitivity; HROM for dual-axis divergence modes; OAPGC for the ecological oral-airway question. Avoid hypermutator/hyper-recombinator causal labels in GTDB and HROM; use “candidate sequence-level” and “candidate architecture-level” language.
2. **Syn2b.** Reserve scaffold/contig-level search, functional-gene-cluster search and structure-based extension plans. The census can be cited as an application enabled by Syn2b anchors, but the primary figure and claims belong to Syn2bANI.
3. **Shared outputs.** Species-level summaries and QC are repository outputs; pair-level GTDB/HROM/OAPGC tables should be deposited on Zenodo with the final manuscript.

## Final figure and table paths

- GTDB joint table: `results/census/gtdb/gtdb_species_ani_sv_joint.tsv`
- GTDB ANI validation: `results/census/gtdb/gtdb_ani_merge_validation.json`
- GTDB quality sensitivity: `results/census/gtdb/effective_repaired/effective_quality_sensitivity.tsv`
- Final GTDB figure: `results/census/gtdb/figures/fig_gtdb_ani_sv_census.pdf`
- Unified structural table: `results/census/unified/species_structural_census_unified.tsv`
