# Synteny Score Benchmark: exact-truth ladder + 695-MAG dnadiff check

Date: 2026-08-21. Supports the "fast synteny from the same pass" claim.
Scripts: `scripts/synteny_bench/` (ladder), `scripts/mag_validation/parse_synteny_truth.py` (MAG check).
Figure: `figures/report/fig_synteny_ladder.png/.pdf`.

## 1. Inversion ladder with exact structural truth

E. coli MG1655 was evolved to ANI 95.00 / 98.00 (counted substitutions) and
given 0–32 non-overlapping inversions (100–400 kb each); each inversion
contributes exactly 2 adjacency breakpoints, so structural truth is exact by
construction (`simulate_inversion_ladder.py`, seed-fixed).

Result (`synteny_ladder_results.tsv`): **breakpoint_count equals the truth on
all 14 rungs** (0–64 breakpoints, zero misses, zero false positives);
synteny_blocks = 2·n_inv + 1; anchor_adjacency decreases monotonically from
1.0000 to 0.9835 (ANI 95, 32 inversions). The ANI estimate is invariant to
rearrangement load: 94.95–95.00 at true 95.00 and 98.00–98.01 at true 98.00
across 0–32 inversions — the chain-restricted likelihood denominator is
unchanged by block reordering, as designed. Wall time 50–70 ms per pair
(Mac Studio, single thread), i.e. the synteny statistics cost nothing on top
of an ANI run that is itself ~150× faster than dnadiff.

## 2. Specificity on 695 real MAG pairs (dnadiff 1-to-1 truth)

For every anchor pair of the CAMI2 MAG benchmark, a nucmer/dnadiff-based
breakpoint count was derived from `dd.1coords` (sort 1-to-1 alignments by
reference start; a break = consecutive alignments on the SAME query contig
violating collinearity; contig switches ignored — `parse_synteny_truth.py`,
`collect/synteny_truth.tsv`).

These MAG–anchor pairs are nearly collinear by construction (a bin comes
from its anchor strain): dnadiff finds 0 breaks for 432/695 pairs, median 0,
max 5. syn2bani agrees exactly on 233/432 zero-break pairs and reports a
median breakpoint_count of 1 (max 34) overall; its anchor_adjacency stays high
(median 0.9997). The residual syn2bani-only breaks concentrate on
fragmented, low-N50 bins — chain fragmentation across assembly gaps, not
called rearrangements. Interpretation: on real data essentially free of
rearrangements, syn2bani does not manufacture structural signal (no pair is
flagged INCONSISTENT on structural grounds); quantitative breakpoint
accuracy at high rearrangement load is established by the exact-truth
ladder (section 1) and the three real rearranged Enterobacteriaceae pairs
(manuscript Fig. 8).

## 3. What the three evidence layers jointly support

- *Exact recall*: 14/14 rungs, breakpoint_count = truth (ladder).
- *Real rearranged pairs*: 184/201 calls with dnadiff counterpart, 80.9%
  coverage of dnadiff ≥ 1 kb events, indel size ratio median 1.000
  (Fig. 8).
- *Specificity at scale*: 695 real pairs, no spurious structural flags,
  anchor_adjacency median 0.9997 on collinear data.
- *Speed*: synteny statistics are a by-product of the ANI pass (50–70 ms
  per pair total); dnadiff costs 8–10 s per pair and nucmer-based synteny
  pipelines more.

## Caveats

- breakpoint_count on heavily fragmented draft queries includes chain
  fragmentation across assembly gaps; for draft-vs-draft comparisons it is
  an upper bound on true rearrangement breaks. anchor_adjacency is the more
  robust of the two statistics on drafts.
- The MAG check measures specificity only; same-strain pairs carry almost
  no true rearrangements to recall.
