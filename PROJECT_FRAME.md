# Syn2bANI project frame

## Objective

Benchmark Syn2bANI as a coordinate-aware ANI and structural-screening method for
complete and draft prokaryotic genomes, and determine what additional biology
the structural layer provides beyond ANI alone.

## Operating modes

| Mode | Input | Primary outputs | Interpretive boundary |
|---|---|---|---|
| M1: complete-genome comparison | FASTA | gated raw ANI, v5-calibrated ANI, AF, structural columns | Calibration was trained only on complete GTDB genomes. |
| M2: draft/MAG comparison | FASTA | raw gated ANI and structural columns | Complete-genome calibration is not used for MAGs. |
| M3: database search | query FASTA versus reference FASTA/sketches | top-hit ANI, structural columns, runtime | skani remains the faster ANI-only screen. |
| M4: catalog census | Syn2b structural layer + independent skani ANI layer | within-species breakpoint and ANI summaries | Not an end-to-end Syn2bANI census. |
| M5: locus-level structural screen | BED or locus coordinates | local chain coverage/SV calls | Requires sufficient local anchors. |

## Benchmark matrix

| ID | Stage | Dataset/scope | Comparator/truth | Required metrics | No-call policy | Status |
|---|---|---|---|---|---|---|
| B01 | exact truth | 12-pair substitution/inversion/indel ladder | skani 0.3.2, FastANI | MAE, bias, flag behavior | `BELOW_DETECTION` is not an estimate | complete |
| B02 | exact-truth robustness | indel, fragmentation and accessory sweeps | self-controls | MAE and AF fidelity | excluded intervals reported | complete |
| B03 | real ANI benchmark | 2,074/2,520-pair v5 training; 43,334-pair held-out GTDB-R207 | ANIm/dnadiff; skani; FastANI; calibrated comparators | MAE, bias, r, genome bootstrap CI | complete-case and declined pairs reported separately | complete |
| B04 | operational mixed regime | 727 high-ANI pairs plus held-out set | raw/calibrated/hybrid; skani; FastANI | complete-case MAE, high-ANI MAE | complete-case denominator explicit | complete |
| B05 | mid-ANI | 15 oral/gut cross-species pairs | ANIm; skani; FastANI | MAE and gate behavior | gate/fallback is the result | complete |
| B06 | MAG | 695 CAMI2 bins | dnadiff ANIm; skani; FastANI | MAE by quality/N50 | raw gated estimate only | complete |
| B07 | database search | 7 queries versus 65,703 GTDB representatives | skani 0.3.2, FastANI, dnadiff follow-up | wall time, RSS, top-hit agreement | FastANI OOM reported; AF floor stated | complete |
| B08 | structural validation | 43,334 held-out pairs; engineered edits | dnadiff/minimap2; constructed controls | correlations, exactness, controls | missed-event modes stated | complete |
| B09 | biological case studies | abfA, FDA-ARGOS E. coli, H. pylori cagPAI | curated status, AMR genotype, lineage stratification | locus AUC/coverage, discordance, associations | genotype, not phenotype; lineage-aware tests | complete |
| B10 | catalog census | HROM, OAPGC, GTDB effective genomes | Syn2b structural + skani ANI layers | coverage, breakpoint fraction, site AUC/enrichment, quality sensitivity | planned/observed and informative denominators reported | complete |

## Comparators and versions

- Syn2bANI 0.1.0 deployed snapshot `fe0f36c`, fragmentation fix `c974f5f`, event-count fix `40b7795`.
- Benchmark-only calibration v6 is not embedded in deployed software.
- skani 0.3.2 for all ANI comparator claims.
- FastANI 1.33/1.34.
- MUMmer/dnadiff 3.23 for ANIm and structural truth.
- Companion Syn2b census/structural snapshot `64717ab`.

## Required provenance

Software commits, model hashes, scripts, input tables, checksums and output
locations are mapped in `REPRODUCIBILITY_MANIFEST.md`. Complete pair-level
census tables require Zenodo deposition before acceptance; aggregate and
verification outputs are in the repository.
