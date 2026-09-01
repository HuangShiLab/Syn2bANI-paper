# What still has to run on the HPC

Ordered by evidential value per CPU-hour. Tasks 1 and 2 need no new alignment and
no new digestion — they are parse passes over files already on disk.

Paths follow `scripts/gtdb50k/run_dnadiff_slice.sh`:

```
WORK=/lustre1/g/aos_shihuang/Syn2bANI-paper/results/gtdb50k
PY=/group/aos_shihuang/conda/bin/python3
ROOT=/lustre1/g/aos_shihuang/Syn2bANI-paper
```

Pull first — the scripts below were added in this commit.

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
```

Then the same over the high-ANI set:

```bash
$PY $ROOT/scripts/compute_dnadiff_events.py \
    --pairs   $WORK/high_ani_pairs_ready.tsv \
    --outdir  $WORK/out_high_ani \
    --outfile $WORK/dnadiff_events_high_ani_all.tsv \
    --workers 32
```

---

## 2. Deduplicate the high-ANI outputs

**Cost: minutes.** `syn2b_inverted_fraction_high_ani_all.tsv` and
`dnadiff_inverted_fraction_high_ani_all.tsv` each carry **195 repeated pairids**;
`high_ani_truth.tsv` carries 65. Merging two tables on a repeated key multiplies
rows, so every count derived from those files without deduplication is inflated —
the >=97% ANIm subset reads 5,655 pairs where 3,826 exist.

`analyze_invfrac_error_model.py` deduplicates defensively at load, so the numbers in
`inverted_fraction_comparison_report.md` are correct. But the cause should be found
rather than papered over — most likely an interrupted run appending to an existing
output instead of rewriting it. Check whether the runner opens its output in append
mode on resume, then regenerate.

```bash
for f in syn2b_inverted_fraction_high_ani_all dnadiff_inverted_fraction_high_ani_all high_ani_truth; do
    echo -n "$f: "
    tail -n +2 $WORK/$f.tsv | cut -f1 | sort | uniq -d | wc -l
done
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
$PY $ROOT/scripts/run_syn2b_inverted_fraction.py \
    --pairs      $WORK/pairs_50k.tsv \
    --genome-dir /lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all \
    --syn2b      /lustre1/g/aos_shihuang/Syn2b/target/release/syn2b \
    --enzymes    BcgI \
    --tgt-dir    $WORK/syn2b_tgts_cache_bcgI \
    --out        $WORK/syn2b_inverted_fraction_50k_bcgI.tsv \
    --workers    4
```

---

## 4. Syn2bANI `breakpoint_count`: fix the reference-side inflation, then re-run

**Cost: a code change, then a re-run of the SV comparison.**

`breakpoint_count = n_chains - n_chained_contigs` subtracts the *query's* contig
term but not the reference's, so a fragmented reference adds exactly `n_ref - 1`
(measured: 10 -> 29 -> 207 junctions at `n_ref` = 1, 20, 200 on an unrearranged
genome). The fix is the positive-contradiction rule already in Syn2b — count a
broken adjacency only when the other genome positively contradicts it, i.e. when an
endpoint already has degree >= 2 there, rather than merely failing to confirm it
(`Syn2b/src/synteny/scoring.rs`). Under that rule a contig break produces no
junction, because a break is an absence of evidence.

Until this lands, `breakpoint_count` correlations against dnadiff carry a
fragmentation term on both sides and `SV_REANALYSIS.md` should be read with that in
mind.

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
