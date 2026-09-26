# HROM hypermode analysis

This analysis adapts the logic of SynTracker Fig. 4 to the completed HROM
within-species census. Instead of SynTracker APSS and inStrain popANI, it uses:

- **skani ANI** as the sequence-similarity axis; and
- the Syn2b `structural` score (conserved landmark adjacencies over their
  union) as the structural-similarity axis.

After filtering pairs with fewer than 250 shared tags or
`observable_fraction < 0.5`, 60,014,738 of 62,872,901 census pairs remained.
The top 5% on each axis contained 3,000,736 pairs. Species enrichment in each
set was tested with one-sided hypergeometric tests and Benjamini-Hochberg FDR;
521 clusters had at least 200 qualifying pairs and were tested.

## Mode counts

| Mode | Species clusters |
|---|---:|
| structural-enriched | 96 |
| ANI-enriched | 160 |
| enriched in both sets | 36 |
| not differentially enriched | 229 |

Overall, skani ANI and the Syn2b structural score were only weakly correlated
(Spearman rho = 0.240 on a fixed 500,000-pair subsample), so the two top-5%
sets emphasize different biological properties.

## Highest-confidence candidate hyper-recombinators

These species are overrepresented in the structural top-5% set relative to the
ANI top-5% set (`log2[structural/ANI] > 0`). They are candidates for elevated
recombination/structural rearrangement.

| Species | total pairs | structural top 5% | ANI top 5% | log2 structural/ANI |
|---|---:|---:|---:|---:|
| *Alloprevotella* sp905369775 | 1,365,109 | 321,121 | 3,835 | 6.39 |
| *Lautropia mirabilis* | 2,382,702 | 328,945 | 3,790 | 6.44 |
| SDRW01 sp007845485 | 363,087 | 94,531 | 5,519 | 4.10 |
| *Prevotella* nanceiensis | 1,727,861 | 221,228 | 10,400 | 4.41 |
| *Actinomyces* graevenitzii | 2,788,689 | 267,212 | 2,059 | 7.02 |
| *Mogibacterium* diversum | 369,179 | 71,034 | 1,392 | 5.67 |
| *Rothia* mucilaginosa | 434,888 | 60,526 | 282 | 7.74 |
| *Corynebacterium* matruchotii | 1,599,189 | 146,079 | 4,152 | 5.14 |

## Candidate sequence-divergent / hypermutator-like taxa

These species are overrepresented in the ANI top-5% set relative to the
structural top-5% set (`log2[structural/ANI] < 0`). Because skani ANI is coarser
than popANI, this should be read as sequence-divergence enrichment, not direct
proof of a hypermutator phenotype.

| Species | total pairs | structural top 5% | ANI top 5% | log2 structural/ANI |
|---|---:|---:|---:|---:|
| *Prevotella* pallens | 2,293,115 | 42,113 | 1,130,017 | -4.75 |
| *Prevotella* pleuritidis | 188,660 | 43,327 | 180,811 | -2.06 |
| F0428 sp003043955 | 186,075 | 21,973 | 155,433 | -2.82 |
| *Peptidiphaga* sp000466165 | 123,753 | 29,547 | 123,086 | -2.06 |
| *Tannerella forsythia* | 54,285 | 7,003 | 52,951 | -2.92 |
| *Porphyromonas gingivalis* | 56,953 | 9,019 | 53,265 | -2.56 |
| *Streptococcus mutans* | 27,261 | 6,270 | 27,138 | -2.11 |
| *Bifidobacterium longum* | 11,781 | 5,258 | 11,476 | -1.13 |

## Files

- `hypermode_species_enrichment.tsv`: full per-cluster statistics.
- `hypermode_figure_points.tsv.gz`: sampled plotting data.
- `hypermode_analysis_report.json`: exact parameters, counts, thresholds and
  correlation.
- `fig_hrom_hypermodes.pdf` / `.png`: SynTracker-style two-panel figure.
