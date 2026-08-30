# SV comparison, re-analysed with confounders controlled

Date: 2026-08-30. Reproduce with `python3 scripts/sv_reanalysis.py results/gtdb50k`
(writes `sv_reanalysis_metrics.tsv`). n = 43,334 pairs.

Supersedes the correlation conclusions in
[SV_COMPARISON_REPORT.md](SV_COMPARISON_REPORT.md) and
[SV_DNADIFF_FILTERED_CORRELATION.md](SV_DNADIFF_FILTERED_CORRELATION.md), which
reported **raw** correlations. A raw correlation cannot separate "both methods
measure rearrangement" from "both methods measure divergence" or "both methods
measure contig count", and at 85–95% ANI on draft assemblies all three are in
play. That is the same failure that made the legacy Syn2b synteny score agree
with SynTracker APSS at r = +0.982 while measuring the wrong thing.

**Contig count is recoverable exactly**, so it can be controlled for. Since
`breakpoint_count = n_chains − n_chained_contigs` and `synteny_blocks = n_chains`:

    n_chained_contigs = synteny_blocks − breakpoint_count

Divergence is controlled with **ANIm from `truth_50k.tsv`**, not Syn2bANI's own
`ani_gated`, which would make the control partly circular.

---

## 1. `breakpoint_count` survives every control — the agreement is real

| pair | raw | \|ANIm | \|n_contigs | **\|both** |
|---|---:|---:|---:|---:|
| `breakpoint_count ~ dnadiff_breakpoints` | 0.465 | 0.503 | 0.472 | **0.534** |
| `breakpoint_count ~ dnadiff_breakpoints_min10000` | 0.489 | 0.531 | 0.493 | **0.553** |
| `breakpoint_count ~ dnadiff_large_indels_min10000` | 0.552 | 0.603 | 0.551 | **0.608** |
| `breakpoint_count ~ mm2_breakpoints` | 0.300 | 0.301 | 0.318 | 0.347 |

The correlation **strengthens** under every control, so it is neither a
divergence artifact nor a fragmentation artifact. `breakpoint_count` correlates
with contig count at only **r = 0.046** (ρ = 0.091): the `n_chains −
n_chained_contigs` subtraction does its job on real draft assemblies. In the
≥95% ANIm subset, `breakpoint_count ~ dnadiff_large_indels_min10000` reaches
**ρ = 0.806**.

This is the column the paper's SV claims should rest on.

## 2. `synteny_blocks` is about half assembly quality — it is not a rearrangement metric

| column | r vs ANIm | r vs n_contigs | ρ vs n_contigs |
|---|---:|---:|---:|
| `breakpoint_count` | 0.197 | **0.046** | 0.091 |
| `synteny_blocks` | 0.305 | **0.767** | 0.652 |
| `mm2_breakpoints` | 0.025 | 0.484 | 0.597 |
| `dnadiff_blocks` | −0.204 | 0.289 | 0.272 |
| `af_query` | **0.449** | 0.106 | 0.070 |
| `anchor_adjacency` | 0.118 | 0.138 | 0.237 |

Contig count averages **31.3** against a `synteny_blocks` mean of **62.2** — so
**50% of "synteny blocks" are just contig starts**, by construction. It correlates
with contig count at r = 0.767, and still at 0.717 inside the ≥95% ANIm subset
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

    dnadiff_breakpoints = 5.35 × breakpoint_count + 290.3      (r = 0.465)

An intercept of **290** is the problem. Among the 3,025 pairs where Syn2bANI
reports **zero** breakpoints, dnadiff still reports a median of **92** (mean 186)
at a median ANIm of 86.5%. Two genomes at 86% ANI do not carry 92 real
rearrangement events.

`dnadiff_blocks` alone explains **R² = 0.600** of `dnadiff_breakpoints`; adding
ANIm and contig count moves it only to 0.605. dnadiff's "Breakpoints" counts
every 1-to-1 alignment boundary, so at this divergence it is dominated by
**alignment fragmentation**, not by biology.

Consequence for how the result is stated: **r ≈ 0.5 is a floor, not a ceiling.**
The agreement is measured against a truth axis that is roughly 60% noise of a
kind Syn2bANI does not share. Reporting it as "moderate agreement" understates
the method; reporting it as validation of a *count* against a clean truth
overstates the reference.

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
(n = 654) is also where the correlations are best (ρ = 0.806), which is
consistent with the diagnosis in §3: the metric looks better where the truth is
cleaner, and that regime is barely sampled.

---

## 5. What to compute next — all from files already on the HPC

**(a) Inverted aligned fraction from the existing `dd.1coords`.** The parser
already exists (it was used for the min5000/min10000 re-parse). Compute

    Σ |E2 − S2| over blocks whose query interval is reversed
    ────────────────────────────────────────────────────────
              Σ |E2 − S2| over all blocks

This is a **length-weighted ratio**, so it has no intercept problem: alignment
fragmentation splits a block into two but leaves the total inverted length
unchanged. It is the direct truth for the upstream Syn2b `inverted_fraction`
signal, which validated against constructed truth at slope 0.968, R² 0.9993 over
a 512× range of event lengths.

**(b) `Inversions`, `Relocations`, `Translocations` from the dnadiff `.report`
files** — dnadiff's own *classified* structural events, far fewer and far closer
to a real event count than raw `Breakpoints`.

**(c) Resample pairs at ≥97% ANI.** Currently n = 2.

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
