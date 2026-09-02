# SV comparison, re-analysed with confounders controlled

Date: 2026-09-02. Reproduce with `python3 scripts/sv_reanalysis.py results/gtdb50k`
(writes `sv_reanalysis_metrics.tsv`). n = 43,334 pairs.

**Correction (2026-09-02):** The `breakpoint_count` implementation previously
subtracted only the query-side contig term (`n_chained_contigs`) and not the
reference-side term. A fragmented reference therefore added `n_ref − 1` spurious
breakpoints. This was fixed in Syn2bANI `c974f5f` and the 43,334 held-out pairs
were re-computed. The numbers below are from the corrected run.

Supersedes the correlation conclusions in
[SV_COMPARISON_REPORT.md](SV_COMPARISON_REPORT.md) and
[SV_DNADIFF_FILTERED_CORRELATION.md](SV_DNADIFF_FILTERED_CORRELATION.md), which
reported **raw** correlations. A raw correlation cannot separate "both methods
measure rearrangement" from "both methods measure divergence" or "both methods
measure contig count", and at 85–95% ANI on draft assemblies all three are in
play.

**Contig count is recoverable exactly**, so it can be controlled for. Since
`breakpoint_count = n_chains − n_chained_contigs` and `synteny_blocks = n_chains`:

    n_chained_contigs = synteny_blocks − breakpoint_count

Divergence is controlled with **ANIm from `truth_50k.tsv`**, not Syn2bANI's own
`ani_gated`, which would make the control partly circular.

---

## 1. `breakpoint_count` after the reference-side fix

| pair | raw | \|ANIm | \|n_contigs | **\|both** |
|---|---:|---:|---:|---:|
| `breakpoint_count ~ dnadiff_breakpoints` | 0.133 | 0.150 | 0.339 | **0.414** |
| `breakpoint_count ~ dnadiff_breakpoints_min10000` | 0.147 | 0.166 | 0.334 | **0.410** |
| `breakpoint_count ~ dnadiff_large_indels_min10000` | 0.255 | 0.281 | 0.380 | **0.453** |
| `breakpoint_count ~ mm2_breakpoints` | 0.064 | 0.061 | 0.309 | 0.340 |

The raw correlation is much lower than before the fix (0.133 vs 0.465) because
the previous value was inflated by the reference-fragmentation term that both
`breakpoint_count` and `dnadiff_breakpoints` shared. After removing it, the
partial correlation controlled for contig count is the relevant signal: **r =
0.339** for `dnadiff_breakpoints` and **0.380** for large indels. Controlling for
both divergence and contig count raises the correlation to **0.414**.

`breakpoint_count` now correlates with query contig count at **r = −0.377**
(ρ = −0.029): the metric is no longer a proxy for assembly fragmentation.

In the ≥95% ANIm subset (n = 654), `breakpoint_count ~
dnadiff_large_indels_min10000` reaches **ρ = 0.674**.

## 2. `synteny_blocks` is about half assembly quality — it is not a rearrangement metric

| column | r vs ANIm | r vs n_contigs | ρ vs n_contigs |
|---|---:|---:|---:|
| `breakpoint_count` | 0.123 | **−0.377** | −0.029 |
| `synteny_blocks` | 0.305 | **0.771** | 0.730 |
| `mm2_breakpoints` | 0.025 | 0.491 | 0.651 |
| `dnadiff_blocks` | −0.204 | 0.262 | 0.297 |
| `af_query` | **0.449** | 0.008 | 0.008 |
| `anchor_adjacency` | 0.118 | 0.005 | −0.072 |

Contig count averages **38.6** against a `synteny_blocks` mean of **62.2** — so
**62% of "synteny blocks" are just contig starts**, by construction. It correlates
with contig count at r = 0.771, and still at 0.759 inside the ≥95% ANIm subset
where divergence cannot explain it.

`SV_DNADIFF_FILTERED_CORRELATION.md` concluded that "breakpoint_count **and
synteny_blocks** are the validated rearrangement metrics … the columns the SV
claims in the paper should rest on." **The `synteny_blocks` half of that
conclusion does not hold**: its agreement with `dnadiff_blocks` (0.494) is in
substantial part two fragmentation measures agreeing with each other, and it
drops to 0.443 once contig count is held constant. It should be reported as an
assembly-structure statistic, not a rearrangement one.

`af_query` is the most divergence-driven column of all (r = 0.449 vs ANIm),
which supports the existing decision to reposition it as coverage/completeness.
`mm2_breakpoints` is a weak third-party check: nearly half its variance tracks
contig count.

## 3. The truth axis is itself contaminated

After the reference-side fix:

    dnadiff_breakpoints = 1.43 × breakpoint_count + 422.0      (r = 0.133)

An intercept of **422** is the problem. Among the 8,312 pairs where Syn2bANI
reports **zero** breakpoints, dnadiff still reports a median of **103** (mean 220)
at a median ANIm of 88.4%. Two genomes at 88% ANI do not carry 103 real
rearrangement events.

`dnadiff_blocks` alone explains **R² = 0.600** of `dnadiff_breakpoints`; adding
ANIm and contig count moves it only to 0.646. dnadiff's "Breakpoints" counts
every 1-to-1 alignment boundary, so at this divergence it is dominated by
**alignment fragmentation**, not by biology.

Consequence for how the result is stated: the low raw correlation is in part
because the reference is noisy. The partial correlation (r = 0.414 after
controlling for ANIm and contig count) is the better estimate of how much true
rearrangement signal `breakpoint_count` captures.

## 4. The validation set is in the wrong ANI range for a strain-level tool

| ANIm | n | share |
|---|---:|---:|
| < 85 | 1,864 | 4.3% |
| 85–90 | 25,747 | 59.4% |
| 90–95 | 15,069 | 34.8% |
| 95–97 | 652 | 1.5% |
| **≥ 97** | **2** | **0.0%** |

This is almost entirely an **inter-species** comparison set. The regime the tool
targets — strain-level, >97% ANI — is represented by two pairs. The ≥95% subset
(n = 654) still shows the strongest signal (`breakpoint_count ~
dnadiff_large_indels_min10000`: ρ = 0.674), consistent with the diagnosis in §3:
the metric looks better where the truth is cleaner, and that regime is barely
sampled.

---

## 5. What to compute next

**(a) Inverted aligned fraction from the existing `dd.1coords` — done.** See
`dnadiff_inverted_fraction.tsv` and `inverted_fraction_comparison_report.md`.
The fixed-reference `raw_inverted_fraction` correlates with dnadiff at r = 0.9355
across the full held-out set.

**(b) `Inversions`, `Relocations`, `Translocations` from the dnadiff `.report`
files — done.** See `dnadiff_events_50k.tsv`.

**(c) Resample pairs at ≥97% ANI.** Currently n = 2. This remains the highest
priority for validating the SV metrics in the strain range the tool is meant for.

### The general principle behind (a)

Every observation process here fragments the genome into segments — contigs
(assembly), 1-to-1 blocks (nucmer), chains (Syn2bANI). Splitting a genome into K
segments removes about K − 1 adjacencies, so **any statistic defined as a count
of transitions picks up a term linear in K**, while a statistic defined as
`Σ(length with property) / Σ(total length)` is invariant to splitting, because a
split preserves both numerator and denominator.

That single fact explains the 290 intercept here, the 50%-contig content of
`synteny_blocks`, and the 119 false junctions a naive adjacency difference
produced on a 120-contig assembly upstream. Where a length-weighted ratio is
available, it should be preferred over a transition count.
