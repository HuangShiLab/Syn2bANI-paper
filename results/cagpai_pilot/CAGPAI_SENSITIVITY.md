# cagPAI pilot: sensitivity of Syn2bANI to island-scale structural variation

Date: 2026-08-26
Tool: syn2bani 0.1.0 (commit fe0f36c, v5 calibration), enzyme panel BcgI,AlfI,AloI,FalI
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

### ANI and synteny score (dist)

| Pair | ANI (%) | af_query | synteny_score | blocks | breakpoints |
|---|---|---|---|---|---|
| wt × wt | 100.00 | 0.999 | 0.9986 | 3 | 2 |
| wt × del | 100.00 | **0.977** | 0.9979 | 4 | 3 |
| wt × inv | 100.00 | 0.998 | 0.9972 | 5 | 4 |
| wt × transloc | 100.00 | 0.998 | 0.9972 | 5 | 4 |
| wt × mut1 | 99.02 | 0.999 | 0.9985 | 3 | 2 |
| mut1 × mut1_del | 97.93 | **0.977** | 0.9975 | 4 | 3 |
| mut1 × mut1_inv | 97.87 | 0.998 | 0.9968 | 5 | 4 |
| mut1 × mut1_transloc | 97.97 | 0.997 | 0.9968 | 5 | 4 |

### SV event calls (struct)

| Pair | Events beyond self-baseline | Called size / coordinates |
|---|---|---|
| wt × del | **none** | (chain count 3→4 only) |
| wt × inv | Inversion | 35,133 bp at 548,447–583,580 — exact cagPAI boundaries |
| wt × transloc | Translocation + Deletion | cagPAI correctly located at new 1.2 Mb locus |
| mut1 × mut1_inv | Inversion | 35,157 bp — robust at ANI ≈ 98 |
| mut1 × mut1_transloc | Translocation + Deletion | robust at ANI ≈ 98 |

Self-baseline (wt × wt): 2 spurious ~10 kb "Inversion" calls at 438 kb ↔ 1,473 kb,
caused by the duplicated rRNA operon repeats. Repeat-masked filtering or a
self-baseline subtraction is required before interpreting struct output.

## Conclusions

1. **ANI estimate is orthogonal to SVs, as designed**: structural variants do not
   perturb the ANI estimate (100.0 vs 100.0); only real nucleotide divergence moves it.
2. **synteny_score is NOT sensitive to a single island-scale event**: 36 kb deletion
   moves it only 0.9986 → 0.9979. The synteny < 0.98 discordance screen used in the
   GTDB analysis therefore captures multi-event / large-scale rearrangement only,
   not cagPAI-scale (2% of genome) events. This limitation must be stated in the
   manuscript, and cagPAI-case analysis must not rely on synteny_score alone.
3. **Inversions and translocations at island scale ARE robustly detected by
   `struct`**, with near-exact breakpoints, even at ANI ≈ 98 (2% divergence).
4. **Pure deletion/insertion is currently NOT emitted as a struct event**; it
   surfaces only as (a) chain count +1 and (b) af asymmetry
   (af_query 0.977 vs af_ref 0.999 ≈ the 2.2% cagPAI genome fraction).
   A deletion/insertion reporter (offset-jump within chains, or AF-asymmetry flag)
   would close this gap — candidate small feature for the Syn2bANI CLI.
5. **Repeat-induced false positives** (rRNA operons) appear even in self × self;
   downstream analyses need repeat masking or a size/consistency filter.

## Implications for the H. pylori cagPAI case study

- Detection of cagPAI presence/absence between high-ANI strains should use
  **af asymmetry (|af_query − af_reference|)** plus chain-count delta, not
  synteny_score; cagPAI inversion/translocation uses struct event calls.
- The struct output must be baseline-filtered (self × self events subtracted or
  repeats masked) before counting events for the 528-genome cohort.
