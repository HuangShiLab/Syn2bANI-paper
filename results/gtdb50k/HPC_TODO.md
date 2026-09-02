# What still has to run on the HPC

Ordered by evidential value per CPU-hour. Tasks 1 and 2 need no new alignment and
no new digestion — they are parse passes over files already on disk.

Paths follow `scripts/gtdb50k/run_dnadiff_slice.sh`:

```
WORK=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k
PY=/group/aos_shihuang/conda/bin/python3
ROOT=/lustre1/g/aos_shihuang/Syn2bANI-paper
```

---

## Status overview

| Task | Status | Notes |
|---|---|---|
| 1. Classified dnadiff events | **done** | `dnadiff_events_50k.tsv` and `dnadiff_events_high_ani_all.tsv` written. |
| 2. Deduplicate high-ANI outputs | **done** | Cause: duplicate rows in `high_ani_pairs_ready.tsv` (pairs sampled into both 95–97 and 97–100 strata), not append mode. Files deduped and stratum corrected. |
| 3. Single-enzyme BcgI pass | **done** | Output `syn2b_inverted_fraction_50k_bcgI.tsv` written (5.7 MB, 43,334 pairs). Validation report: `BCGI_ERROR_MODEL_VALIDATION.md`. |
| 4. `breakpoint_count` reference-side inflation | **done** | Fix pushed to Syn2bANI `main` (`c974f5f`). HPC binary rebuilt; wrapper `3974244` completed; `s2b_50k.tsv` regenerated. Reports updated. |
| 5. Closed-genome control | **not started** | Low priority; kept for reviewer response. |

---

## 1. Classified dnadiff events — the count channel is still unvalidated

**Cost: minutes. No alignment, no digestion.** `dd.report` already exists for every
pair; only its `[Feature Estimates]` block has never been read.

Everything validated so far is the *orientation* channel, and it is validated
against a single scalar. dnadiff separates events by kind — `Inversions`,
`Relocations`, `Translocations`, and `Insertions` (which are indels, not
rearrangements). Two questions turn on that split:

- Does Syn2b's junction count track *rearrangements*, or does it track indels? The
  intercept in `dnadiff_breakpoints = 5.35 * breakpoint_count + 290.3`
  (`SV_REANALYSIS.md`) is unexplained. If it is mostly `Insertions`, the truth axis
  is contaminated and the whole breakpoint comparison needs restating.
- Detection power was measured in simulation (`Syn2b/docs/PHASE2_DETECTION_POWER.md`)
  as ~2.6 kb L50 for inversions and ~1.5 kb for translocations on BcgI. Classified
  counts let that be checked against real genomes rather than simulated ones.

```bash
$PY $ROOT/scripts/compute_dnadiff_events.py \
    --pairs   $WORK/pairs_50k.tsv \
    --outdir  $WORK/out \
    --outfile $WORK/dnadiff_events_50k.tsv \
    --workers 32

$PY $ROOT/scripts/compute_dnadiff_events.py \
    --pairs   $WORK/high_ani_pairs_ready.tsv \
    --outdir  $WORK/out_high_ani \
    --outfile $WORK/dnadiff_events_high_ani_all.tsv \
    --workers 32
```

**Done.** 50k: 43,334 pairs, all complete. High-ANI: 7,710 pairs, 7,710 complete
(after deduplication; see Task 2).

---

## 2. Deduplicate the high-ANI outputs

**Cost: minutes.** Done.

The original `high_ani_pairs_ready.tsv` contained 65 duplicate `pairid`s because
the same pairs were sampled into both the 95–97 and 97–100 strata. This propagated
into `syn2b_inverted_fraction_high_ani_all.tsv`,
`dnadiff_inverted_fraction_high_ani_all.tsv`, `high_ani_truth.tsv`, and
`dnadiff_events_high_ani_all.tsv`. The runners do **not** open outputs in append
mode; the duplication was in the input pair list.

Fix applied on HPC:

- Deduplicated `high_ani_pairs_ready.tsv`, keeping the stratum that matches the
  measured ANIm (`<97` → 95–97, `>=97` → 97–100).
- Re-added the `pairid` column to the cleaned pair list.
- Deduplicated the four output tables by `pairid` and filtered to the cleaned pair
  list.

Current duplicate counts:

```bash
for f in syn2b_inverted_fraction_high_ani_all dnadiff_inverted_fraction_high_ani_all high_ani_truth dnadiff_events_high_ani_all high_ani_pairs_ready; do
    echo -n "$f: "
    tail -n +2 $WORK/$f.tsv | cut -f1 | sort | uniq -d | wc -l
done
# all report 0
```

---

## 3. Single-enzyme BcgI, same pairs — a sharp falsifiable prediction

**Cost: one full digestion + comparison pass over the held-out set.**

Everything reported so far used the four-enzyme panel `BcgI,AlfI,AloI,FalI`
(`run_syn2b_inverted_fraction.py`), but a real 2bRAD library is usually one enzyme.
Since the error is now a closed-form function of the shared-landmark count `m`,

```
SE = sqrt( 1.504 * p(1-p) / m + 0.0205^2 )
```

dropping to one enzyme should cut `m` by roughly 4x and raise the SE by exactly the
amount that formula predicts. If it does, the error model is right and the panel
choice can be made on paper for any future enzyme set. If the SE rises by more, the
model is missing a term — most likely because a sparser panel also raises the
detection floor, which would be worth knowing before the paper claims a design rule.

`--enzymes` is new in this commit; **give it its own `--tgt-dir`**, since the cache
keys on accession alone and would otherwise silently reuse four-enzyme TGTs.

```bash
# Submitted as scripts/gtdb50k/s12_syn2b_bcgI_invfrac.slurm
# Job: 3974133, partition amd, 16 CPUs, 4 h.
$PY $ROOT/scripts/run_syn2b_inverted_fraction.py \
    --pairs      $WORK/pairs_50k.tsv \
    --genome-dir /lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all \
    --syn2b      /lustre1/g/aos_shihuang/Syn2b/target/release/syn2b \
    --enzymes    BcgI \
    --tgt-dir    $WORK/syn2b_tgts_cache_bcgI \
    --out        $WORK/syn2b_inverted_fraction_50k_bcgI.tsv \
    --workers    16
```

---

## 4. Syn2bANI `breakpoint_count`: fix the reference-side inflation, then re-run

**Cost: a code change, then a re-run of the SV comparison.**

`breakpoint_count = n_chains - n_chained_contigs` subtracted the *query's* contig
term but not the reference's, so a fragmented reference added exactly `n_ref - 1`
(measured: 10 → 29 → 207 junctions at `n_ref` = 1, 20, 200 on an unrearranged
genome).

**Fix applied in `Syn2bANI/src/core/chain_ani.rs` (`c974f5f`).** A transition
between two chains along the query is now counted only when the reference genome
positively contradicts the query adjacency — at least one endpoint has two
neighbours along its reference contig. A reference contig break therefore produces
no breakpoint, because a break is an absence of evidence. Unit tests added for
reference-fragmentation invariance and inversion counting.

**Re-run in progress.** The fixed release binary was copied to HPC and
`s2b_out/*.tsv` was cleared so the new `breakpoint_count` is written for all
43,334 held-out pairs. The original 190-task array hit the account
`MaxSubmitJobsPerAccount` limit, so a single-job wrapper
(`scripts/gtdb50k/s14_s2b_single_wrapper.slurm`, job `3974197`) is running the
190 slices locally in parallel with 8 workers. After it finishes, regenerate:

```bash
$PY $ROOT/scripts/gtdb50k/analyze_sv_comparison.py   # or whatever merges s2b_out
```

and update `SV_REANALYSIS.md` / `SV_COMPARISON_REPORT.md` with the corrected
correlations.

---

## 5. Closed-genome control — low priority, and here is why

The obvious way to separate the `sigma0 = 0.0205` floor from assembly fragmentation
is to restrict to pairs where both genomes are single-contig. That is worth doing
eventually, but the evidence already points away from fragmentation as the cause:
after removing the sampling term, the residual correlates with `observable_fraction`
at **+0.004**, and controlling for it moves the overall correlation from 0.9355 to
0.9354. The floor also *falls* below 0.0205 at high ANI (SD 0.0122 at >=99.5% ANIm),
which is a divergence effect, not a contig-count effect. Run this only if a reviewer
asks.

---

## Not needed any more

**Resampling pairs at >=97% ANIm.** This was on the list because held_out_50k has
n = 2 there. The high-ANI set already covers it: 3,826 pairs at >=97% measured ANIm,
with slope 1.0063, r = 0.9960 and SD 0.0135. Combined with held_out_50k's 80–97%
coverage, the ANIm axis is continuous from 80% to 100% with no gap. See the banded
tables in `inverted_fraction_comparison_report.md`.
