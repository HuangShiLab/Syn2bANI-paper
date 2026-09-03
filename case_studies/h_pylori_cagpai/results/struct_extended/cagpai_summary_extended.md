# H. pylori cagPAI extended state analysis

## Overview

We classified 528 H. pylori genomes from Song et al. (2026) for cagPAI state using two complementary Syn2bANI outputs:

1. **Presence/absence of 28 marker CDS** (HP0520–HP0547 from strain 26695) → `empty`, `partial`, `complete`.
2. **`syn2bani struct --bed` of each genome vs 26695** → structural variation inside the cagPAI coordinates (NC_000915.1: 547,327–583,481).

Combining both gives four extended states:

- `empty` — 0–few markers present (island absent).
- `partial` — some but not all markers present (island degraded/deleted).
- `complete_collinear` — all markers present and no inversion/translocation inside cagPAI.
- `complete_rearranged` — all markers present but at least one inversion or translocation inside cagPAI.

## Results

### Overall distribution (n = 528)

| State | Count | Percentage |
|---|---:|---:|
| empty | 85 | 16.10% |
| partial | 11 | 2.08% |
| complete_collinear | 12 | 2.27% |
| complete_rearranged | 420 | 79.55% |

The striking result is that **most strains that retain the cagPAI carry a detectable rearrangement** relative to the 26695 reference. Only ~2% are fully collinear.

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

The extended state distribution differs significantly across Correa’s cascade stages (χ² = 27.91, df = 9, p = 9.9 × 10⁻⁴):

| Group | n | empty | partial | complete_collinear | complete_rearranged |
|---|---:|---:|---:|---:|---:|
| NAG | 167 | 39 (23.35%) | 1 (0.60%) | 3 (1.80%) | 124 (74.25%) |
| AG | 133 | 14 (10.53%) | 3 (2.26%) | 1 (0.75%) | 115 (86.47%) |
| IM | 110 | 23 (20.91%) | 1 (0.91%) | 3 (2.73%) | 83 (75.45%) |
| GC | 118 | 9 (7.63%) | 6 (5.08%) | 5 (4.24%) | 98 (83.05%) |

GC shows the highest proportion of both `complete_collinear` and `partial` states, while NAG has the highest `empty` rate.

### Association with FastBAPS lineage

FastBAPS lineage is strongly associated with cagPAI state (χ² = 58.68, df = 12, p = 3.9 × 10⁻⁸):

| Lineage | n | empty | partial | complete_collinear | complete_rearranged |
|---|---:|---:|---:|---:|---:|
| L2 | 197 | 8 (4.06%) | 7 (3.55%) | 2 (1.02%) | 180 (91.37%) |
| L3 | 149 | 41 (27.52%) | 1 (0.67%) | 2 (1.34%) | 105 (70.47%) |
| L4 | 48 | 3 (6.25%) | 1 (2.08%) | 3 (6.25%) | 41 (85.42%) |
| L5 | 113 | 26 (23.01%) | 2 (1.77%) | 5 (4.42%) | 80 (70.80%) |
| L6 | 21 | 7 (33.33%) | 0 (0.00%) | 0 (0.00%) | 14 (66.67%) |

L2 is dominated by rearranged-but-present islands, whereas L3, L5 and L6 have much higher empty rates.

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

- `syn2bani struct --bed query.fna hp26695.fna -o query.vs_hp26695.bed`
- cagPAI region on NC_000915.1: **547,327–583,481** (from cag1 start to cagA end; ±2 kb buffer used for overlap).
- For strains already called `complete` by markers:
  - any INV/TRA overlapping the region → `complete_rearranged`
  - any DEL ≥10 kb overlapping the region → `partial`
  - otherwise → `complete_collinear`
- Strains called `partial` or `empty` by markers keep that state.

## Files

- `cagpai_states_extended.tsv` — per-genome extended state and list of cagPAI-overlapping SVs.
- `cagpai_association.tsv` — contingency tables and χ² tests.
- `figures/cagpai_state_by_*.png` — stacked-bar visualizations.

## Interpretation for the paper

1. **Syn2bANI adds information that ANI alone cannot provide.** Two strains can have nearly identical ANI but differ in cagPAI collinearity; only the structural call reveals this.
2. **H. pylori cagPAI is highly dynamic.** The reference-collinear form is rare in this global collection; rearranged forms are the norm.
3. **Disease-stage differences are detectable.** GC-enriched `complete_collinear` and `partial` states suggest that both island presence and island architecture may matter for gastric-cancer risk.
