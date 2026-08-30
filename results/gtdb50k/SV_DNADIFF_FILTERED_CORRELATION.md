# dnadiff gap-filtered SV truth vs Syn2bANI metrics (43,334 GTDB held-out pairs)

> **Correction (2026-08-30).** The correlations below are **raw**, and conclusion
> 2 does not survive controlling for assembly fragmentation. Contig count is
> recoverable as `synteny_blocks − breakpoint_count`; on that axis
> `synteny_blocks` correlates at **r = 0.767** (0.717 even within the ≥95% ANIm
> subset) while `breakpoint_count` correlates at only **0.046**. Half of
> `synteny_blocks` is contig starts by construction. `breakpoint_count` remains
> validated — its agreement *strengthens* under every control — but
> **`synteny_blocks` should not be presented as a rearrangement metric.** See
> [SV_REANALYSIS.md](SV_REANALYSIS.md).

Date: 2026-08-26
Inputs: `sv_truth_50k_min5000.tsv`, `sv_truth_50k_min10000.tsv` (re-parsed
dd.1coords with min-gap 5k/10k; HPC job dd_filter 3945832),
`s2b_50k.tsv` (syn2bani dist --verbose).

## Correlations (n = 43,334)

| dnadiff truth | Syn2bANI metric | min 5k: pearson / spearman | min 10k: pearson / spearman |
|---|---|---|---|
| breakpoints | anchor_adjacency | −0.04 / −0.26 | −0.05 / −0.27 |
| breakpoints | breakpoint_count | +0.48 / +0.54 | +0.49 / +0.56 |
| breakpoints | synteny_blocks | +0.50 / +0.61 | +0.49 / +0.61 |
| blocks | breakpoint_count | +0.44 / +0.60 | +0.44 / +0.60 |
| blocks | synteny_blocks | +0.49 / +0.62 | +0.49 / +0.62 |
| large_indels | breakpoint_count | +0.51 / +0.66 | **+0.55 / +0.70** |
| large_indels | anchor_adjacency | +0.00 / −0.17 | −0.04 / −0.27 |
| large_indels | af_query | +0.08 / +0.11 | +0.04 / +0.09 |

## Conclusions

1. **Gap-filtering dnadiff does not rescue anchor_adjacency.** anchor_adjacency is
   AF-like (fraction of genome covered by chains) and is insensitive to
   rearrangement count — consistent with the cagPAI pilot (36 kb deletion
   moved it only 0.9986 → 0.9979). It should be repositioned in the
   manuscript as a genome-coverage / completeness metric, not a
   rearrangement metric.
2. **breakpoint_count and synteny_blocks are the validated rearrangement
   metrics**: spearman 0.54–0.70 against dnadiff breakpoints/blocks/
   large-indels, improving with the 10k filter (large-scale events).
   These are the columns the SV claims in the paper should rest on.
3. Residual decorrelation is expected: dnadiff (nucmer) and Syn2bANI
   (restriction-tag anchors) differ in resolution and in draft-assembly
   contig handling; minimap2 orthogonal truth (s11, pending) will provide
   a third point of comparison.
