# cagPAI extended-state re-analysis after circular-origin filtering

**Date:** 2026-09-06

## What changed

The original `syn2bani struct` pipeline reported a `complete_rearranged` state for
420 / 528 (79.6%) *H. pylori* cohort genomes versus reference 26695. A diagnostic
re-analysis showed that 400 / 420 (95.2%) of these carried a structural call
spanning >50% of the ~1.67 Mb chromosome, most commonly `TRA:1459-1666206` (283
genomes, 99.8% of the chromosome).

Such genome-spanning calls are not biological rearrangements; they arise because
bacterial chromosomes are circular but assemblies start at arbitrary coordinates.
When the query assembly uses a different arbitrary start than the reference, the
linear SV detector splits the collinear chromosome into two chains and reports a
translocation covering the whole replicon. Because the cagPAI window
(547,327–583,481) lies inside this span, the overlap test was always true and
inflated `complete_rearranged`.

## Fix

Syn2bANI `struct` now supports:

- `--circular <contig>[,...]` — declare reference/query contigs as circular.
- `--artifact-threshold <fraction>` — filter calls whose reference or query span
  exceeds this fraction of the contig length (default 0.5).

The H. pylori runner now calls:

```bash
syn2bani struct --bed --circular NC_000915.1 <query> hp26695.fna -o <out>.bed
```

Verified on GCA_000521245.1: the `TRA:1459-1666206` call is removed while local
SVs are retained.

## New extended-state counts

| State | n | % |
|---|---:|---:|
| empty | 85 | 16.10 |
| partial | 11 | 2.08 |
| complete_collinear | 145 | 27.46 |
| complete_rearranged | 287 | 54.36 |

`complete_rearranged` dropped from 420 to 287 after removing circular-origin
artifacts. The presence/absence axis (`empty`/`partial`/`complete_*`) is
unaffected by the filter.

## Association with metadata

### Marginal associations (not lineage-adjusted)

| Metadata field | χ² | df | p |
|---|---:|---:|---:|
| Disease stage (group) | 24.60 | 9 | 0.0034 |
| FastBAPS lineage | 58.75 | 12 | 3.8 × 10⁻⁸ |

### Crude and lineage-stratified 2×2 associations

The table below gives crude odds ratios (Woolf 95% CI) and Cochran–Mantel–Haenszel
statistics stratified by FastBAPS lineage. Only the first contrast (cagPAI
presence, GC vs NAG) had a nominally significant crude association; the other
three were not significant before stratification.

| Contrast | n | crude OR (95% CI) | crude χ² p | CMH p | OR_MH | Breslow-Day p |
|---|---:|---:|---:|---:|---:|---:|
| cagPAI presence: GC vs NAG | 285 | 2.16 (1.13–4.13) | 0.018 | 0.121 | 1.70 | 0.005 |
| cagPAI rearrangement: GC vs NAG | 230 | 1.11 (0.64–1.93) | 0.700 | 0.574 | 1.17 | 0.153 |
| cagPAI presence: GC/IM vs AG/NAG | 528 | 1.14 (0.73–1.78) | 0.576 | 0.172 | 1.39 | 0.021 |
| cagPAI rearrangement: GC/IM vs AG/NAG | 432 | 1.02 (0.68–1.52) | 0.928 | 0.951 | 1.01 | 0.780 |

For the presence GC-vs-NAG contrast, the crude OR of 2.16 (p = 0.018) is
attenuated to 1.70 after lineage stratification (CMH p = 0.12). The Breslow-Day
test is significant (p = 0.005), indicating that the odds ratio is not uniform
across FastBAPS strata. Cell counts in the key discordant cells are small:
non-cancer cagPAI-negative genomes total only 15 across the cohort, and after
stratification they are distributed as 4 / 9 / 1 / 1 / 0 across FastBAPS L2–L6.
The association is therefore compatible with a modest, lineage-heterogeneous
effect that this design is under-powered to resolve.

### Artifact-threshold sensitivity

The rearrangement association was recomputed under artifact-span thresholds from
0.20 to 0.50. The number of `complete_rearranged` genomes ranged from 254 to 287,
but the crude and CMH p-values for rearrangement vs disease stage remained
non-significant across the entire range (crude p = 0.51–0.70 for GC vs NAG;
Fig. S12a). The presence association is unchanged because presence is determined
by marker alignment, not by structural calls.

### Reverse-complement diagnostic

A whole-genome reverse-complement assembly convention would appear as a single
large inversion spanning most of the chromosome. Among 497 genomes with an
unfiltered call spanning ≥20% of the chromosome, the largest such call was a
translocation in 496 cases and an inversion in only 1 case. The dominant pattern
is therefore circular-origin *rotation* (different arbitrary start coordinate),
not strand-direction convention. This supports the use of `--circular` together
with a span threshold rather than an inversion-specific filter.

## Consequences for the manuscript

- The cagPAI case study cannot claim a direct disease-stage association for
  rearrangement state; that signal was largely a circular-origin artifact and
  does not survive filtering or lineage stratification.
- The presence/absence axis (`empty`/`partial`/`complete`) remains associated
  with lineage and geography and is still biologically meaningful; it should be
  the focus if the case study is retained.
- The crude presence association (GC vs NAG) is attenuated, not abolished, after
  lineage stratification. It should be reported with both crude and CMH
  statistics, the Breslow-Day heterogeneity test, and the small discordant-cell
  counts, rather than as a simple "no longer significant" statement.
- The structural-comparison story should pivot to the GTDB-R207 evidence
  (`raw_inverted_fraction`, `breakpoint_count` with contig control) and the
  engineered cagPAI pilot panel, where the exact 36,154 bp deletion is recovered.
