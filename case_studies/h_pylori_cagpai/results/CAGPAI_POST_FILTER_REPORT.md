# cagPAI extended-state re-analysis after circular-origin filtering

**Date:** 2026-09-04

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

### Lineage-stratified associations (Cochran–Mantel–Haenszel by FastBAPS)

| Contrast | CMH χ² | df | p | OR_MH |
|---|---:|---:|---:|---:|
| cagPAI presence: GC vs NAG | 1.40 | 1 | 0.236 | 1.60 |
| cagPAI rearrangement: GC vs NAG | 0.80 | 1 | 0.370 | 1.29 |
| cagPAI presence: GC/IM vs AG/NAG | 1.55 | 1 | 0.213 | 1.39 |
| cagPAI rearrangement: GC/IM vs AG/NAG | 0.47 | 1 | 0.495 | 1.15 |

After stratifying by FastBAPS lineage, neither cagPAI presence nor cagPAI
rearrangement state is significantly associated with disease stage. The marginal
disease-stage association was therefore driven by lineage structure, which is
itself confounded with geography and cohort composition.

## Consequences for the manuscript

- The cagPAI case study can no longer claim a direct disease-stage association
  for rearrangement state.
- The presence/absence axis (`empty`/`partial`/`complete`) remains associated
  with lineage and geography and is still biologically meaningful; it should be
  the focus if the case study is retained.
- The structural-comparison story should pivot to the GTDB-R207 evidence
  (`raw_inverted_fraction`, `breakpoint_count` with contig control) and the
  engineered cagPAI pilot panel, where the exact 36,154 bp deletion is recovered.
