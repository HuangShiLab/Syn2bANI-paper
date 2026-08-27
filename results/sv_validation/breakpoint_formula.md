# Breakpoint-formula validation (`synteny_stats`, `src/core/chain_ani.rs`)

Date: 2026-08-14. Code: Syn2bANI main @ 98177dc + the fix described below
(working tree, not committed). Binary: `target/release/syn2bani` rebuilt
after the fix. All runs: default 4-enzyme panel (BcgI,AlfI,AloI,FalI),
mismatch tolerance 2, `--verbose`. Raw tables: `{simindel,simacc,simfrag}{,_fixed}.tsv` in this directory.

## What the formula counts

`synteny_stats(chains, anchors)` computes, per pair:

- `possible` = within-query-contig anchor adjacencies, Σ_c (k_c − 1) over
  contigs with k_c anchors;
- `conserved` = adjacencies that are consecutive inside one chain,
  Σ (chain_len − 1) = (chained anchors) − n_chains;
- `breakpoint_count = possible − conserved`.

**Before the fix**, `possible` ranged over *all* seed anchors. Since greedy
chain extraction rejects anchors that do not lie on a high-scoring collinear
path (multi-mapping repeats, off-diagonal mismatch-tolerance hits, runs
shorter than `min_chain_anchors = 4`), every rejected anchor inflated
`breakpoint_count` by one. The count therefore decomposed as

```
breakpoint_count = (n_chains − n_chained_contigs)   # genuine chain transitions
                 +  n_unchained_anchors             # chaining-rejected anchors
```

and the second term dominated by two orders of magnitude: an E. coli pair at
95% ANI with **zero** rearrangements reported ~670 "breakpoints"
(~15% of anchors are rejected by chaining).

**After the fix**, `possible` ranges over chained anchors only, so
`breakpoint_count = n_chains − n_chained_contigs` exactly — the number of
chain-to-chain transitions along the query. A clean inversion (+/−/+ chain
structure on one contig) gives exactly 2. `anchor_adjacency` uses the same
chained-only denominator and is now ~1.0 for collinear pairs.

The INCONSISTENT flag is untouched: `unreliable()` still receives the old
all-anchor statistic, now named `SyntenyStats::unconserved`
(`possible_all − conserved`), because the 0.5-per-anchor threshold in
`results/gating_flag/RULES.md` was calibrated on exactly that quantity
(AUC 0.803 — its power comes largely from the rejected-anchor term).
`flag`, `ani*`, `n_chains`, and every other output column are bit-identical
before/after on all 38 validation rows below (asserted by script).

## Empirical validation on constructions with known rearrangements

### simindel — one clean 400 kb inversion per genome (12 genomes, ANI 0.85–0.999)

Expected: exactly 2 breakpoints (the inversion's two junctions).

| genome | true ANI | bp before | bp after | chains | verdict |
|---|---|---|---|---|---|
| q_ani0.8500 | 0.85 | 148 | 3 | 4 | +1: sparse anchors (retention 0.10) split one collinear arm |
| q_ani0.8800 | 0.88 | 281 | **2** | 3 | exact |
| q_ani0.9000 | 0.90 | 232 | **2** | 3 | exact |
| q_ani0.9200 | 0.92 | 564 | **2** | 3 | exact |
| q_ani0.9400 | 0.94 | 613 | **2** | 3 | exact |
| q_ani0.9500 | 0.95 | 710 | **2** | 3 | exact |
| q_ani0.9600 | 0.96 | 728 | **2** | 3 | exact |
| q_ani0.9700 | 0.97 | 832 | **2** | 3 | exact |
| q_ani0.9800 | 0.98 | 835 | **2** | 3 | exact |
| q_ani0.9900 | 0.99 | 797 | 3 | 4 | +1 chain split |
| q_ani0.9950 | 0.995 | 737 | 5 | 6 | +3 chain splits |
| q_ani0.9990 | 0.999 | 774 | **2** | 3 | exact |

The +/−/+ 3-chain structure of the inversion is recovered correctly at every
level; excess counts come from extra chains, i.e. the adaptive `max_skip`
chain-break test splitting collinear arms at anchor-sparse regions — a
chain-continuity artifact, not a miscount of structure.

### simacc — 5 accessory (non-homologous, length-preserving) blocks, core ANI 95%

Expected: 0 for acc0.00; 5 for the rest (each block interrupts collinearity
once — one chain transition per block).

| genome | accessory frac | bp before | bp after | chains | verdict |
|---|---|---|---|---|---|
| acc0.00 | 0.0 | 671 | **0** | 1 | exact |
| acc0.10 | 0.1 | 692 | **5** | 6 | exact |
| acc0.20 | 0.2 | 676 | **5** | 6 | exact |
| acc0.30 | 0.3 | 539 | **5** | 6 | exact |
| acc0.40 | 0.4 | 701 | **5** | 6 | exact |
| acc0.50 | 0.5 | 613 | **5** | 6 | exact |

### simfrag — draft fragmentation (20/50/100/200 contigs, ±50% reverse-complemented, shuffled order), ANI 95%

Expected: 0 — fragmentation and contig orientation create no within-contig
rearrangement, and contig boundaries are not counted by construction
(each contig's first chain is free). `flipped` vs `fwd` rows isolate the
orientation effect.

| genome | contigs | bp before | bp after | chains | verdict |
|---|---|---|---|---|---|
| q95_c20 | 20 | 715 | 2 | 22 | 2 contigs split into 2 chains |
| q95_c20_fwd | 20 | 714 | 2 | 22 | identical with flips disabled |
| q95_c50 | 49 | 716 | 2 | 51 | same |
| q95_c50_fwd | 49 | 714 | 2 | 51 | same |
| q95_c100 | 99 | 717 | 1 | 99 | 1 split |
| q95_c100_fwd | 99 | 713 | 1 | 99 | same |
| q95_c200 | 201 | 692 | **0** | 200 | exact |
| q95_c200_fwd | 201 | 689 | **0** | 200 | exact |

Reverse-complementation adds **zero** breakpoints (flipped vs fwd identical
to within ±1), confirming strand-canonical matching works; the residual 0–2
is contigs whose tag-poor halves split into two chains.

## Verdict

**The formula was miscounting — bug, fixed.** The old `breakpoint_count`
conflated genuine chain transitions with chaining-rejected anchors
(100–800 artifact counts on clean E. coli pairs, scaling with anchor count,
not with rearrangement). The structural term underneath was always correct.

Fix (working tree, uncommitted):

- `src/core/chain_ani.rs` — `synteny_stats` now counts `possible`
  adjacencies over chained anchors only
  (`breakpoint_count = n_chains − n_chained_contigs`); the old all-anchor
  difference is kept as `SyntenyStats::unconserved` and is still what
  `unreliable()` thresholds, preserving the calibrated INCONSISTENT flag
  exactly. `anchor_adjacency` denominator likewise chained-only.
  Doc comments updated (`SyntenyStats`, `FLAG_MAX_BP_PER_ANCHOR`,
  `unreliable`, `ChainAniResult::{breakpoint_count,anchor_adjacency,unreliable}`,
  `src/cli/ani.rs` header).
- Unit tests: strengthened `synteny_stats_single_inversion` (exact count),
  added `synteny_stats_inversion_two_breakpoints` (+/−/+ = 2),
  `synteny_stats_ignores_unchained_anchors`,
  `synteny_stats_fragmented_query`. `cargo test --lib`: **72 passed**
  (69 pre-existing + 3 new), green in both debug and release.
- Drive-by fix, pre-existing: debug-mode `attempt to multiply with overflow`
  in the `seqs()` test helper (plain `*` on an LCG constant) made 7 lib
  tests panic under plain `cargo test` at 98177dc too; changed to
  `wrapping_mul` (release behaviour unchanged, sequences identical).

Residual, characterized not fixed (correct-by-design behaviour):
`breakpoint_count` counts *chain transitions*, so anything that splits a
chain — including anchor-sparse regions at low retention and deletions
larger than `max_gap` — adds to the count. On clean data this is 0–3
spurious transitions; the classical SV signal (2 per inversion, 1 per
accessory block) sits cleanly on top. Equal-length homologous replacements
that keep tags phased do not break chains and are not counted (documented
limitation), as are chains shorter than 4 anchors.
