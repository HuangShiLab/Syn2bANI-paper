# H. pylori cagPAI extended state analysis

## Overview

We classified 528 H. pylori genomes from Song et al. (2026) for cagPAI state using two complementary Syn2bANI outputs:

1. **Presence/absence of 28 marker CDS** (HP0520–HP0547 from strain 26695) → `empty`, `partial`, `complete`.
2. **`syn2bani struct --bed` of each genome vs 26695** → structural variation inside the cagPAI coordinates (NC_000915.1: 547,327–583,481).

Combining both gives four extended states:

- `empty` — 0–few markers present (island absent).
- `partial` — some but not all markers present (island degraded/deleted).
- `complete_collinear` — all markers present and no local inversion/translocation inside cagPAI.
- `complete_rearranged` — all markers present but at least one local inversion or translocation inside cagPAI.

Because *H. pylori* has a circular chromosome, a different arbitrary choice of
start coordinate between two assemblies can produce a single genome-spanning SV
that overlaps every locus by construction. Such calls were flagged and excluded
from the rearrangement count (genome-spanning = >50% of the 1,667,825 bp
chromosome).

## Results

### Overall distribution after circular-origin filtering (n = 528)

| State | Count | Percentage |
|---|---:|---:|
| empty | 85 | 16.10% |
| partial | 11 | 2.08% |
| complete_collinear | 145 | 27.46% |
| complete_rearranged | 287 | 54.36% |

Before filtering, 133 additional complete-marker genomes were classified as
`complete_rearranged` solely because of a genome-spanning circular-origin call.
After filtering, **about half of the genomes that retain cagPAI carry a local
rearrangement** relative to the 26695 reference, while the other half are
reference-collinear once coordinate-system differences are removed.

### Validation with engineered controls

We ran the same pipeline on the nine control strains used in the cagPAI pilot study:

| Control | Presence/absence | Extended state | Interpretation |
|---|---|---|---|
| wt | complete | complete_collinear | intact, collinear |
| hp26695 | complete | complete_collinear | intact, collinear |
| mut1 | complete | complete_collinear | intact, collinear |
| inv | complete | complete_rearranged | inversion inside cagPAI |
| transloc | complete | complete_rearranged | translocation inside cagPAI |
| mut1_inv | complete | complete_rearranged | inversion inside cagPAI |
| mut1_transloc | complete | complete_rearranged | translocation inside cagPAI |
| del | empty | empty | island deleted |
| mut1_del | empty | empty | island deleted |

All nine controls are classified consistently with their known biology.

### Association with disease stage

The unfiltered extended state distribution differs significantly across Correa’s
cascade stages (χ² = 27.91, df = 9, p = 9.9 × 10⁻⁴), but this association is
confounded by FastBAPS lineage. After circular-origin filtering and
lineage stratification (Cochran–Mantel–Haenszel), the disease-stage association
is no longer significant:

| Contrast | CMH χ² | df | p | OR_MH |
|---|---:|---:|---:|---:|
| cagPAI presence (GC vs NAG) | 1.4019 | 1 | 0.2364 | 1.599 |
| cagPAI rearrangement (GC vs NAG) | 0.8028 | 1 | 0.3703 | 1.294 |
| cagPAI presence (advanced vs early) | 1.5539 | 1 | 0.2126 | 1.389 |
| cagPAI rearrangement (advanced vs early) | 0.4664 | 1 | 0.4946 | 1.152 |

Filtered disease-stage counts:

| Group | n | empty | partial | complete_collinear | complete_rearranged |
|---|---:|---:|---:|---:|---:|
| NAG | 167 | 39 (23.35%) | 1 (0.60%) | 45 (26.95%) | 82 (49.10%) |
| AG | 133 | 14 (10.53%) | 3 (2.26%) | 37 (27.82%) | 79 (59.40%) |
| IM | 110 | 23 (20.91%) | 1 (0.91%) | 29 (26.36%) | 57 (51.82%) |
| GC | 118 | 9 (7.63%) | 6 (5.08%) | 34 (28.81%) | 69 (58.47%) |

### Association with FastBAPS lineage

FastBAPS lineage remains strongly associated with cagPAI state after filtering
(χ² = 58.75, df = 12, p = 3.8 × 10⁻⁸):

| Lineage | n | empty | partial | complete_collinear | complete_rearranged |
|---|---:|---:|---:|---:|---:|
| L2 | 197 | 8 (4.06%) | 7 (3.55%) | 66 (33.50%) | 116 (58.88%) |
| L3 | 149 | 41 (27.52%) | 1 (0.67%) | 30 (20.13%) | 77 (51.68%) |
| L4 | 48 | 3 (6.25%) | 1 (2.08%) | 21 (43.75%) | 23 (47.92%) |
| L5 | 113 | 26 (23.01%) | 2 (1.77%) | 26 (23.01%) | 59 (52.21%) |
| L6 | 21 | 7 (33.33%) | 0 (0.00%) | 2 (9.52%) | 12 (57.14%) |

L2 has the lowest `empty` rate and the highest `complete_collinear` rate,
whereas L3, L5 and L6 have much higher `empty` rates. The lineage signal is
robust to circular-origin correction.

### Association with country and phylogenetic population

Country-level and population-level patterns mirror the FastBAPS result (e.g., high empty rates in Sweden and hspNEurope; high rearranged rates in China and hspEAsia). Some cells are small, so chi-square is not reported, but the stacked-bar visualizations are provided.

## Methods

### Marker-based classification

- Reference: H. pylori 26695 (NC_000915.1).
- Markers: 28 CDS loci from HP0520 (cag1) to HP0547 (cagA).
- Each marker was queried against the 528 assemblies with `minimap2 -x asm5`.
- A marker was considered present if an alignment covered ≥80% of its length at ≥80% identity.
- State thresholds: `empty` = ≤10% markers present; `partial` = 10–90%; `complete` = ≥90%.

### Structural classification

- `syn2bani struct --bed --circular NC_000915.1 query.fna hp26695.fna -o query.vs_hp26695.bed`
- cagPAI region on NC_000915.1: **547,327–583,481** (from cag1 start to cagA end; ±2 kb buffer used for overlap).
- For strains already called `complete` by markers:
  - any INV/TRA overlapping the region → `complete_rearranged`
  - any DEL ≥10 kb overlapping the region → `partial`
  - otherwise → `complete_collinear`
- Strains called `partial` or `empty` by markers keep that state.
- Complete-marker strains whose only cagPAI-overlapping SV spans >50% of the
  chromosome are reclassified as `complete_collinear` (circular-origin artifact).

## Files

- `cagpai_states_extended.tsv` — per-genome extended state and list of cagPAI-overlapping SVs.
- `cagpai_association.tsv` — contingency tables and χ² tests.
- `figures/cagpai_state_by_*.png` — stacked-bar visualizations.

## Interpretation for the paper

1. **Syn2bANI adds information that ANI alone cannot provide.** Two strains can have nearly identical ANI but differ in cagPAI presence/absence or local collinearity; only the structural call reveals this.
2. **H. pylori cagPAI is lineage-structured, not disease-stage associated.** After correcting for circular-origin artifacts and stratifying by FastBAPS lineage, the disease-stage association disappears. The robust signal is lineage-specific presence/absence (`empty` rates 4–33% across lineages).
3. **Circular-origin filtering is essential for circular chromosomes.** Without it, a coordinate-system difference between two assemblies is misread as a genome-spanning rearrangement that overlaps every locus, inflating apparent rearrangement rates.
