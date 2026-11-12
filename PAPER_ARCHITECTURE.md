# Syn2bANI paper architecture

## Reader logic

The manuscript moves from method design to controlled accuracy, real benchmark
accuracy, practical search, metric-specific structural validation, biological
loci, and catalog-scale population screens.

## Section-to-claim map

| ID | Section | One-sentence claim | Primary display | Evidence |
|---|---|---|---|---|
| A1 | Syn2bANI estimates ANI and architecture | Fixed Type IIB anchors support chain-restricted ANI and coordinate-aware structural outputs in one pass | Fig. 1 | Methods; algorithm figure |
| A2 | Exact truth | Syn2bANI is accurate under controlled substitutions and robust to indels, fragmentation and accessory content | Fig. 2 | `results/simulations/` |
| A3 | GTDB benchmark | Calibrated Syn2bANI outperforms shipped comparators on a genome-disjoint ANIm benchmark, but calibration boundaries are explicit | Fig. 3; Table 1 | `results/gtdb50k/` |
| A4 | Search/runtime | Syn2bANI is slower than skani but much faster than dnadiff for SV follow-up and returns ANI plus structure | ED Fig. 1; Table 2 | `results/one_to_all/` |
| A5 | Structural validation | Structural metrics are metric-specific: inversion fraction is strongly validated; breakpoint count is moderate; adjacency/chains are diagnostics | Fig. 4 | `results/gtdb50k/`; `results/sv_validation/` |
| A6 | Functional loci | Coordinate-carrying anchors reveal high-ANI strains with abfA deletion or resistance-genotype/structural differences | Fig. 5 | abfA and FDA-ARGOS case studies |
| A7 | cagPAI | Syn2bANI structural calls resolve cagPAI architecture, but associations require lineage-aware interpretation | Fig. 6 | H. pylori case study |
| A8 | HROM census | HROM shows non-monotonic ANI-architecture relationships and candidate dual-axis phenotypes | ED Fig. 3 | `results/census/hrom/` |
| A9 | OAPGC census | Cross-site conspecifics show ANI-derived distance enrichment, not a global breakpoint excess | ED Fig. 4 | `results/census/oapgc/` |
| A10 | GTDB census | The approach scales to hundreds of millions of effective pairs and shows declining breakpoint prevalence at very high ANI | ED Fig. 5 | `results/census/gtdb/` |

## Display policy

Main text: Figures 1–6 and Tables 1–2. Extended Data: Figures 1–5 and Tables
1–2. Secondary sensitivity analyses, detailed cohort tables and historical
process records are in the supplement or repository, not main-text prose.
