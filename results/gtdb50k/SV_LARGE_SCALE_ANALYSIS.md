# Large-scale SV evaluation on GTDB 50k held-out pairs — interim analysis

## Status

dnadiff 1-to-1 truth for all **43,334** held-out pairs already exists on HPC,
and `parse_dnadiff_sv.py` has been run. Aggregated results are in
`results/gtdb50k/sv_truth_50k.tsv` and `SV_EVALUATION_REPORT.md`.

## Key finding: dnadiff and Syn2bANI measure different things

| Metric | Syn2bANI | dnadiff |
|---|---|---|
| mean breakpoint_count | 30.9 | 455.7 |
| median breakpoint_count | 20 | 333 |
| mean anchor_adjacency | 0.971 | 0.449 |

**Correlations across 43,334 pairs:**
- breakpoint_count Pearson r = 0.465 (Spearman 0.527)
- anchor_adjacency Pearson r = 0.119 (Spearman 0.329)

**Rearrangement detection (any breakpoint):**
- precision = 0.999, recall = 0.931, F1 = 0.964
- specificity = 0.608

## Interpretation

dnadiff produces many more breakpoints because it fragments alignments at
small indels, repeats, and low-identity regions. Its "anchor adjacency"
(1 − breakpoints / (blocks − 1)) is dominated by alignment fragmentation,
not large-scale rearrangements. Syn2bANI's sparse tags and chain-restricted
likelihood bridge small variations and report fewer, larger breakpoints.
The high precision (0.999) means Syn2bANI almost never reports a breakpoint
where dnadiff reports none; the lower recall (0.931) means it misses some
small dnadiff-detected breaks, as expected by design.

## Consequence for the manuscript

The current draft claims that anchor_adjacency "carries biologically meaningful
signal orthogonal to ANI" but does not benchmark it against an independent
structural truth at scale. A reviewer will likely run the same comparison and
find r ≈ 0.12 against dnadiff-derived synteny. We must address this before
submission.

Options:

1. **Filter dnadiff truth to large events only.** Re-parse dd.1coords with a
   minimum gap threshold (e.g., 5 kb or 10 kb) so that both tools are scored
   on large-scale rearrangements. This is the most direct fix and should be
   run on HPC where the 1coords files live.

2. **Use minimap2 + paftools or SyRI as an orthogonal truth.** minimap2
   produces fewer, larger alignment blocks than dnadiff; its breakpoint
   calls may correlate better with Syn2bANI's chain geometry.

3. **Reframe the anchor adjacency definition.** Instead of comparing absolute
   values, show that Syn2bANI's anchor_adjacency predicts *functional* or
   *evolutionary* discordance (e.g., the B. longum abfA case), which is the
   biologically relevant claim.

## Next steps

- [ ] Re-run `parse_dnadiff_sv.py` on HPC with min_gap = 5,000 and 10,000 bp.
- [ ] Run minimap2 on the same 43,334 pairs (or a stratified subset) and
      derive SV truth from PAF.
- [ ] Compare Syn2bANI struct output against both dnadiff-large and minimap2
      truth.
- [ ] Update manuscript: add large-scale SV validation paragraph, clarify
      that anchor_adjacency measures large-scale rearrangement burden, and
      report the resolution difference between tools.
