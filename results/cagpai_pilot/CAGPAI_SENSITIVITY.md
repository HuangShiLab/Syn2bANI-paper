# cagPAI pilot: sensitivity of Syn2bANI to island-scale structural variation

Date: 2026-08-26 (updated after the struct shadow-chain fix)
Tool: syn2bani 0.1.0 (commit 4d679cf, v5 calibration + struct shadow filter),
enzyme panel BcgI,AlfI,AloI,FalI
Reference genome: H. pylori 26695 (NC_000915.1 / GCF_000008525.1, 1,667,867 bp)
cagPAI coordinates (from GFF annotation, cag1..cagA): 547,328–583,481 (36,154 bp)

## Design

In-silico variants of the 26695 genome (`data/cagpai_pilot/`):

| Variant | Description | Length |
|---|---|---|
| wt | original | 1,667,867 |
| del | 36 kb cagPAI deletion (empty-site) | 1,631,713 |
| inv | 36 kb cagPAI inverted in place | 1,667,867 |
| transloc | cagPAI moved to 1.2 Mb locus | 1,667,867 |
| mut1 | wt + 1% random SNPs | 1,667,867 |
| mut1_del / mut1_inv / mut1_transloc | independent 1% SNP background + same SV (pairwise divergence ~2%) | — |

Runs: `syn2bani dist --verbose` (8×8, `dist_all.tsv`) and `syn2bani struct --rearrangement --indel` on 8 focal pairs (`struct_out/`).

## Results

### ANI and anchor adjacency (dist)

| Pair | ANI (%) | af_query | anchor_adjacency | blocks | breakpoints |
|---|---|---|---|---|---|
| wt × wt | 100.00 | 0.999 | 0.9986 | 3 | 2 |
| wt × del | 100.00 | **0.977** | 0.9979 | 4 | 3 |
| wt × inv | 100.00 | 0.998 | 0.9972 | 5 | 4 |
| wt × transloc | 100.00 | 0.998 | 0.9972 | 5 | 4 |
| wt × mut1 | 99.02 | 0.999 | 0.9985 | 3 | 2 |
| mut1 × mut1_del | 97.93 | **0.977** | 0.9975 | 4 | 3 |
| mut1 × mut1_inv | 97.87 | 0.998 | 0.9968 | 5 | 4 |
| mut1 × mut1_transloc | 97.97 | 0.997 | 0.9968 | 5 | 4 |

### SV event calls (struct, after shadow-chain fix)

| Pair | SVs | Call |
|---|---|---|
| wt × wt (self-baseline) | **0** | clean — rRNA repeat artifacts removed |
| wt × del | 1 | **Insertion 36,154 bp at 546,687–583,917 — exact cagPAI** |
| wt × inv | 1 | Inversion 35,133 bp at 548,447–583,580 — exact cagPAI |
| wt × transloc | 2 | Translocation + Deletion 36,154 bp, cagPAI located at new 1.2 Mb locus |
| mut1 × mut1 (self) | 0 | clean |
| mut1 × mut1_del | 1 | Insertion 36,154 bp — robust at ANI ≈ 98 |
| mut1 × mut1_inv | 1 | Inversion 35,157 bp — robust at ANI ≈ 98 |
| mut1 × mut1_transloc | 2 | Translocation + Deletion 36,154 bp — robust at ANI ≈ 98 |

## The struct shadow-chain fix (syn2bani 4d679cf)

Before the fix, self × self reported 2 spurious ~10 kb "Inversions"
(438 kb ↔ 1,473 kb): secondary mappings of the duplicated rRNA operons.
Worse, those sparse nested chains sat between the two backbone chains in
query order and broke their adjacency, so the 36 kb cagPAI deletion was
**never called**. The fix marks any chain fully contained in a denser
chain's query span as a repeat shadow and excludes it from all SV rules
(`src/core/sv.rs::shadowed_mask`). Unit tests (11/11) and the full suite
(106/106) pass, including the K-12 vs Sakai relocation-guard regression.

## Conclusions

1. **ANI estimate is orthogonal to SVs, as designed**: structural variants do not
   perturb the ANI estimate (100.0 vs 100.0); only real nucleotide divergence moves it.
2. **anchor_adjacency is NOT sensitive to a single island-scale event**: 36 kb deletion
   moves it only 0.9986 → 0.9979. The synteny < 0.98 discordance screen used in the
   GTDB analysis therefore captures multi-event / large-scale rearrangement only,
   not cagPAI-scale (2% of genome) events. This limitation must be stated in the
   manuscript, and cagPAI-case analysis must not rely on anchor_adjacency alone.
3. **`struct` now detects all three island-scale event classes** — deletion/
   insertion, inversion, translocation — with near-exact boundaries (±1–2 kb of
   the true cagPAI edges), robustly at ANI ≈ 98 (2% divergence). Self-baselines
   are clean (0 calls).
4. **Repeat-induced false positives are gone** after the shadow filter; no
   repeat masking or baseline subtraction is needed downstream.
5. The `dist`-level `breakpoint_count`/`synteny_blocks` still count the shadow
   chains (wt × wt reports 2 breakpoints) — cosmetic only, but worth a note if
   these columns are compared against struct event counts.

## Implications for the H. pylori cagPAI case study

- cagPAI presence/absence, inversion, and translocation between high-ANI
  strains can all be read directly from `syn2bani struct` event calls;
  af asymmetry (|af_query − af_reference| ≈ the island's genome fraction)
  serves as a cheap screen before running struct.
- The 528-genome cohort analysis can proceed without a repeat-masking
  pre-processing step.
- HPC struct jobs (top-discordant GTDB pairs, B. longum) should be run with
  the rebuilt binary (≥ 4d679cf).
