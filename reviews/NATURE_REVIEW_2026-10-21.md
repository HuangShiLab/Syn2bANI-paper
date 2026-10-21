# Nature Methods review round, 21 October 2026

## Review process

Three isolated reviewer runs received the same packet
(`reviews/NATURE_REVIEW_PACKET_2026-10-21.md`) and different emphasis briefs.
The frozen reports are `2026-10-21_R1.md`, `2026-10-21_R2.md` and
`2026-10-21_R3.md`. This file is the post-review synthesis and revision record.

## Consensus strengths

- The restriction-anchor concept is distinctive and useful for joint ANI and
  structural reporting.
- Exact-truth simulations, genome-disjoint GTDB calibration, negative controls,
  reliability flags and quality-sensitivity analyses form a strong evidence base.
- The manuscript now separates deployed v5 from benchmark-only v6 and reports
  uncomfortable comparator results explicitly.
- The three catalog analyses are valuable scale demonstrations when interpreted
  as associative screens.

## Blocking concerns and revision response

### R1-M1. Hybrid threshold selected on the reported benchmark

The reviewer correctly noted that the 98% hybrid threshold was selected on the
same unified set used for the headline score. The manuscript now reports the
unified hybrid result as an operational-threshold sensitivity analysis on
40,629 complete-case pairs. The confirmatory claim is limited to the genome-
disjoint raw and v6-calibrated results. Supplementary Note 6 and
`results/gtdb50k/CALIBRATION_PROVENANCE.md` now state the selection and
evaluation overlap explicitly.

### R2-M1. OAPGC ANI distance is not genome-wide substitution load

The estimand is now described as alignable-sequence ANI-derived distance. The
text states that skani conditions on shared sequence and that contiguity,
aligned fraction, accessory content, ancestry and public-label error can affect
the comparison. The conclusion is now a public-label-associated difference,
not a habitat process, migration direction or SNP-load claim. A renamed
machine-readable table, `species_effects_all_pairs_ani_distance.tsv`, replaces
the historical `SNP_*` terminology without breaking the original file.

### R2-M2. Breakpoint enrichment is not a complete structural census

The structural-validation section now calls breakpoint count a detection
statistic with moderate correlation to dnadiff and states known missed-event
classes. HROM labels are now candidate anchor-statistic enrichment phenotypes
rather than divergence modes. The phrase "at least two detected breakpoints" is
used where the detection limitation matters.

### R2-M3. GTDB cluster scope, missing breakpoints and dominance

New `results/revision/gtdb_cluster_scope_sensitivity.tsv` compares complete
species, all clusters with at least 100 expected pairs, partial clusters at
coverage thresholds, analyses excluding *E. coli*, unweighted medians,
pair-weighted fractions and genus/family aggregation. The manuscript now
distinguishes the 1,311 complete-species median (61.3%) from the pair-weighted
fraction (56.2%) and from all 1,316 eligible clusters (61.6% unweighted; 77.9%
pair-weighted). GTDB estimates are explicitly descriptive observational
summaries without hierarchical uncertainty.

### R3-M1. Exact software and artifact provenance

New `REPRODUCIBILITY_MANIFEST.md` records Syn2bANI commits for the deployed
0.1.0 snapshot, fragmentation fix, 0.1.1 event-count fix and HPC benchmark
snapshot; companion Syn2b commit `64717ab`; v5/v6 model checksums; result-to-
script mapping; and checksums for major analysis artifacts.
It also distinguishes the provided runtime specification from a byte-for-byte
lockfile or container.

### R3-M4. Review-stage data availability

Data availability now separates the immediately available reduced review
verification set from the complete pair-level Zenodo archive planned before
acceptance. It adds the new AF, cluster-scope and ANI-distance files, references
the checksummed manifest, and states the controlled-access route for restricted
collaborator metadata.

### R3-M5. Structural validation was not uniform across metrics

The section and Figure 5 title now use metric-specific wording. Strong agreement
is claimed for the companion-tool inverted fraction; breakpoint count is
moderate and exact only for constructed edits; aligned fraction is a coverage
descriptor rather than SV validation. FDA-ARGOS classification thresholds,
strand/overlap rules, cs-string indel verification and artifact handling are now
documented in Supplementary Table S7.

## Other major changes

- Comparator calibration is now described as separate post-hoc recalibration,
  not identical retraining of Syn2bANI v6; shipped, calibrated and
  FastANI-feature tiers are distinguished.
- The benchmark target is defined as skani-screenable representative same-genus
  pairs, with explicit complete-case and all-pair limitations.
- New AF-stratified benchmark metrics cover strict, low-AF and very-low-AF
  truth tiers (`results/revision/af_stratified_benchmark.tsv`).
- Uncertainty scope is narrowed: genome bootstrap intervals are reported for
  primary Syn2bANI-vs-skani contrasts, not every band or phylum.
- HROM all-pair and high-quality labels receive equal analytical status, and
  examples are restricted to high-quality-retained labels.
- Efficiency definitions were reconciled: 484 ordered comparisons include self;
  3.7/2.8 s rows are ANI plus structural columns; separate struct adds 4.36 s;
  5,000-genome screen/refinement is 406.7 s.

## Minor changes

- Units were changed to ANI points in the abstract and key benchmark text.
- Figure and table mappings were corrected for the unified high-ANI set and
  Enterobacteriaceae structural set.
- Table 1 now gives complete-case denominators and explains the blank high-ANI
  unified cells.
- Table 2 now explains ordered self comparisons and corrects memory values.
- The *S. aureus* country caveat now uses 351/2,211 pairs.
- Supplementary Figure S4 now states that extreme fragmentation increases error.
- The eight-genome cagPAI panel now explicitly includes the matched mut1
  control.
- Figure 8 separates primary disease-stage tests from the exploratory locality
  analysis.
- FDA-ARGOS now states that genotype, not phenotype, was measured.
- The historical 50k shorthand remains only where filenames require it, with
  exact counts in the text.
