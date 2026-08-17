# Short-contig rescue — validation report (2026-08-17)

Status: IMPLEMENTED and deployed (code repo commit 0aabd0c, file
`src/core/chain_ani.rs` only, +571/−4). Fixes the fragmentation-driven upward
bias documented in V8_MLE_VALIDATION.md §3.7 and the main report §2.4.

## 1. Mechanism diagnosis (quantified)

The upward drift is a **miss-count selection**, not a chain-acceptance
problem per se. On the Sakai→K-12 ladder (regenerated with
`prototype/fragment.py` into `prototype/ladder_sakai/`, true ANI constant
98.0535):

- At N50 5 kb, 529 of 1,103 contigs (2.44 Mb, 2,366 of 7,064 tags) carried no
  accepted chain and were dropped from the likelihood; 256 of them (1.15 Mb)
  had ≥2 anchors and were plainly homologous.
- Decomposing the partial estimators: `ani_from_hist` is flat across the
  ladder (98.78 intact → 98.78 at N50 5 kb) while `ani_from_loss` collapses
  (97.82 → 98.76). Diverged tags don't anchor → anchor-poor regions fail
  `min_chain_anchors = 4` or fall outside chains' anchor-bounded fill spans →
  the miss count erodes faster than the hit count → the gamma fit has no
  heterogeneity left to correct and reads high.
- Simply lowering `min_chain_anchors` to 2 barely helped (98.707 → 98.697):
  fill only counts tags *between* a chain's outer anchors; the misses live in
  the tails and unchainable contigs.
- A per-tag join against the intact run (by tag sequence) proved two distinct
  tail populations: tail tags the intact run *also* counts (87% misses — real
  divergence that must enter the likelihood) and tail tags it excludes via
  pass-2 `max_skip` chain breaks (hyper-divergent runs whose inclusion
  overshoots −0.84 at N50 20 kb). Any correct fix must separate these two.

## 2. Design

A third stage inside `fill()`, strictly gated so long-contig assemblies take
a bit-identical path:

- A query contig with `< 8 × min_chain_anchors` in-panel tags (default 32 ≈
  25 kb with the default panel; the smallest 200 kb-rung contig still has 153
  tags) and **no accepted chain** is rescued if it has a *collinear group of
  ≥ 2 unique anchors* (`uniq` = tag occurs once in each genome; the group must
  contain ≥ 1 anchor with ≤ 1 mismatch — a solely-2-mismatch basis can be a
  paralog whose true locus was occurrence-filtered, measured on a fragmented
  self-comparison). Ties between equal-size groups break on total mismatch
  count, not hash-map order (a real rc-control regression).
- Counting replicates the pass-2 chain rule on a complete genome: tags
  **between** the outer basis anchors are counted (hits → histogram, no match
  → miss); a **tail tag** past an outer anchor is counted only while a
  bracketing anchor (from any contig — typically the flanking one) follows
  within `(max_skip+1) × mean tag spacing` bp along the reference. That is
  exactly "would a chain have spanned this locus on the intact genome".
- Rescued contigs contribute AF spans with the same half-median-gap extension
  rule as chains; `synteny_score`/`breakpoint_count`/`n_chains`/SV `chains`
  output are untouched (rescue creates no chains).
- `ChainAniConfig::short_contig_rescue` (default true) disables it for A/B.

## 3. Validation

**Sakai ladder** (`ani`, drift vs intact 98.0535):

| rung | before | after | drift before | drift after | AF before → after |
|---|---|---|---|---|---|
| complete | 98.0535 | 98.0535 | — | — **bit-identical** | 0.7688 → 0.7688 |
| N50 500k | 98.0863 | 98.0863 | +0.03 | +0.03 **identical** | 0.7646 → 0.7646 |
| N50 200k | 98.0808 | 98.0808 | identical | identical | 0.7590 → 0.7590 |
| N50 100k | 98.1239 | 98.1239 | identical | identical | 0.7499 → 0.7499 |
| N50 50k | 98.1947 | 98.1941 | +0.141 | +0.141 | 0.734 → 0.734 |
| N50 20k | 98.3364 | 98.0218 | +0.283 | **−0.032** | 0.683 → 0.724 |
| N50 10k | 98.5440 | 98.0338 | +0.491 | **−0.017** | 0.620 → 0.690 |
| N50 5k | 98.7071 | 98.2544 | +0.654 | **+0.200** | 0.434 → 0.580 |

Cross-tool on the same N50 5k file: skani 97.84 (drift −0.24; complete
98.08), FastANI 97.41/97.13 (drift −0.64/−0.92 depending on direction).
syn2bani's |0.20| now beats both. (On this regenerated ladder skani/FastANI
drift *downward*; the +0.35–0.38 figures in V8 §3.7 were from an earlier
ladder build.)

**draftbench** (8 real ENA *E. coli* drafts): all 8 rc self-controls =
99.9999 before and after. Vs K-12:

| assembly | contigs | before | after | skani | fastani |
|---|---|---|---|---|---|
| GCA_001075925 (contaminated) | 8,025 | 97.9250 | **97.3922** | 97.01 | 97.02 |
| GCA_001077875 | 143 | 98.5083 | 98.5066 | 98.48 | 98.43 |
| GCA_001283205 | 225 | 98.5969 | 98.5927 | 98.28 | 98.30 |
| GCA_001283245 | 121 | 98.5673 | 98.5664 | 98.31 | 98.34 |
| GCA_001283605 | 162 | 98.5601 | 98.5394 | 98.32 | 98.27 |
| GCA_001284145 | 137 | 97.1014 | 97.0743 | 98.41 | 98.34 |
| GCA_001283865 | 88 | 98.5413 | 98.5373 | 96.48 | 96.69 |
| GCA_001284645 | 95 | 97.1401 | 97.1237 | 96.98 | 96.60 |

Worst case (8,025 contigs): bias vs skani **+0.91 → +0.36**. The
GCA_001284145/GCA_001283865 disagreements with skani are pre-existing and
unaffected.

**Regression**: simindel default-panel output **bit-identical** (`diff`
empty; rescue never triggers on long-contig genomes). `cargo test --release`:
101 green (90 lib + 11 integration), including 5 new tests:
`nearest_anchor_distance_both_directions`,
`rescue_basis_requires_two_unique_collinear_anchors`,
`rescue_basis_prefers_exact_group_over_paralog_tie`,
`fragmented_genome_does_not_inflate_estimate` (digest-level: rescue counts
strictly more tags, ends closer to truth than rescue-off, agrees with the
intact estimate), `fragmented_identical_genomes_stay_at_one` (rescue must not
manufacture misses/mismatches; AF only ever increases).

## 4. Caveats

- **Recalibration advised for `--calibrate`/v4**: the v4 linear model's
  features (`af_query`, `af_reference`, `retention`, `n_anchors`, `n_chains`,
  `n_tags_in_chains`) shift on fragmented inputs only: at N50 5 kb `af_query`
  0.434→0.580, `n_tags_in_chains` 4,457→5,523, `retention` 0.9246→0.9106;
  `n_anchors`/`n_chains` unchanged. Complete genomes are bit-identical
  including calibrated output (97.7976 both). If the deployed model is used
  on drafts, retrain it on post-fix features.
- `std_err` shrinks slightly on drafts (more counted tags) — correctly so.
- Residual +0.20 at N50 5 kb comes from contigs with 0–1 placeable anchors
  (~1.3 Mb), which carry no positional evidence; rescuing them would require
  guessing homology and would risk the downward noise the bracket /
  unique-anchor rules were built to exclude.
- Extreme drafts whose contigs *all* fail to chain still return NaN (the
  pass-1 `chains.is_empty()` early return is unchanged) — the rescue only
  runs when at least one real chain exists. Not hit by any benchmark here.
- One benign exception to "synteny untouched": at N50 20 kb pass-2 picked 231
  vs 233 blocks because the pass-1 fit now sees rescued tags, slightly
  shifting `max_skip`.
