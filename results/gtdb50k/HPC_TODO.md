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
| 3. Single-enzyme BcgI pass | **running** | SLURM job `3974133` (`g50k_bcgI_invfrac`) on `amd`, 16 CPUs, 4 h. |
| 4. `breakpoint_count` reference-side inflation | **code done, re-run in progress** | Fix pushed to Syn2bANI `main` (`c974f5f`). New binary copied to HPC; s2b array resubmitted to regenerate `s2b_50k.tsv`. |
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

## 6. Which taxa actually carry large inversions — annotate what we already have

**Cost: minutes.** No new alignment, no new digestion, no downloads.

The inverted-fraction channel is only interpretable on near-closed assemblies.
Contigs are deposited in arbitrary orientation, so a fragmented pair drifts to
`min(f, 1-f) ~ 0.5` carrying no biology. Measured on the >=97% ANIm pairs, with
contig count recovered from `observable_fraction ~ 1 - (K-1)/S`:

| K (estimated) | n | median min(f,1-f) | fraction < 0.05 |
|---|---|---|---|
| 1-2 (closed) | 371 | 0.2880 | 0.332 |
| 3-5 | 250 | 0.3776 | 0.100 |
| 6-10 | 257 | 0.3991 | 0.035 |
| 11-25 | 691 | 0.4281 | 0.023 |
| 26-100 | 1255 | 0.4548 | 0.007 |
| >100 | 1002 | 0.4726 | 0.004 |

The 0.5 plateau is contig-orientation randomness and it disappears exactly as the
assemblies close up. So the 371 near-closed pairs are the only ones in hand where
the number means something — 96 are collinear (< 0.01) and 183 sit at 0.30-0.51.

What is missing is which species they are; only phylum is in the pair tables. This
joins GTDB taxonomy and ranks species by how often their near-closed pairs disagree
in orientation, so the target organism for the inversion case study is chosen from
data rather than guessed.

```bash
$PY $ROOT/scripts/annotate_closed_inversion_pairs.py \
    --results  $WORK \
    --metadata /lustre1/g/aos_shihuang/data/gtdb-r207/metadata/bac120_metadata_r207.tsv \
    --out      $WORK/closed_inversion_pairs.tsv
```

It also re-checks `K_est` against the metadata's real `contig_count` and
`ncbi_assembly_level`, which is the honest version of the table above. Send back
`closed_inversion_pairs.tsv` and the printed species ranking.

---

## 7. The SynTracker isolate cohorts — a controlled test with published answers

**Run this before task 7b.** Task 7b picks organisms from GTDB on structural
grounds; this one uses cohorts whose answer is already published.

Enav, Paz & Ley, *Nat Biotechnol* 43:773-783 (2024),
doi:10.1038/s41587-024-02276-2, pairs SynTracker (synteny) with inStrain (SNPs)
on four isolate collections and reports which mode dominates in each. That gives
a design our own data cannot: **cohorts where the expected answer is already
published, including a positive and a negative control.**

| Cohort | n | Published result | What our `inverted_fraction` must do |
|---|---|---|---|
| *Streptomyces rimosus* M527, different fermentations | 20 isolates, 185 pairs | popANI spans only 0.99990-1.0 (clonal) while APSS spans ~0.90-1.0 — variation is **structural, not SNP** | **positive control**: must show signal where SNP-based tools see nothing |
| hypermutator *E. coli*, 4 mice colonised with 2 ancestral substrains | 185 pairs | inStrain calls *none* the same strain; SynTracker calls *all* the same strain — variation is **SNP, not structural** | **negative control**: must stay near 0 |
| *Neisseria gonorrhoeae* clinical isolates, antibiotic resistant | 12 isolates, 66 pairs | Spearman rho = 0.985 between the two scores — **both modes** | both channels should move together |
| *Helicobacter pylori* clinical isolates, 6 participants | 77 isolates, 21-91 pairs each | mixed; participants 322, 326, 439 carry subpopulations that only one tool calls same-strain | the interesting case, and it connects to the cagPAI cohort |

A negative control is the part we have never had. Everything measured so far shows
the estimator tracks dnadiff; nothing shows it stays quiet when there is nothing to
find. The hypermutator *E. coli* set is exactly that test, and it is a hard one —
those genomes carry a heavy SNP load, which is what strips landmarks.

### Reference genomes (taken from the paper's Methods, not recalled)

| Species | Accession |
|---|---|
| *H. pylori* | GenBank `CP032479.1` |
| *S. rimosus* M527 | `GCF_000331185.2` (ASM33118v2) |
| *N. gonorrhoeae* | `GCF_900087635.2` |

### Getting the reads

The per-isolate SRA accessions are in the paper's Supplementary Tables 2-5, which
are **not in the main PDF** — download the Supplementary Information from the
article page first. Source studies, if the SI is awkward to parse:

- *H. pylori* — Wilkinson, Dickins, Robinson & Winter, *Gut Microbes* 14, 2152306 (2022)
- hypermutator *E. coli* — Ramiro, Durão, Bank & Gordo, *PLoS Biol* 18, e3000617 (2020)
- gut metagenomes (the 1,133-individual analysis) — Suzuki et al., *Science* 377, 1328-1332 (2022)

`scripts/syntracker_validation/02_download_reads.sh` and `03_assemble_array.sh`
already do read download and SPAdes assembly, so the pipeline exists.

### The methodological catch, and how the paper solves it

These are **SPAdes assemblies from short reads**, and task 6 shows the
inverted-fraction channel is uninterpretable on fragmented input. The paper hits
the same wall and works around it:

> the contigs of each assembly were ordered using the reference genomes described
> above and then aligned against each other

Mauve's contig mover orients every contig against the reference, which removes the
arbitrary-orientation artifact. **We need the same preprocessing step**, and it has
a consequence worth stating in the paper rather than discovering in review:
reference-guided orientation biases each contig toward collinearity with the
reference, so inversions whose breakpoints fall on contig boundaries are absorbed
rather than detected. The bias is toward the null — it costs sensitivity, not
specificity — which is acceptable for a positive/negative control design but must
be declared.

Run each cohort both ways (raw contigs, and reference-ordered) and report both. The
gap between them is itself a measurement of how much the orientation artifact costs.

### A claim of ours that this cohort decides

`MATH_REVIEW.md` section 6 records our measurement that APSS responds to
divergence rather than order (delta -0.00016 structure / -0.08751 divergence),
and concludes agreement with APSS is not a goal. The *S. rimosus* cohort is a
direct test of that: popANI is pinned at ~1.0 there, so if APSS tracked only
divergence it would be flat, and the paper reports it spanning 0.90-1.0. Either
our measurement generalises less far than stated, or something else separates the
two settings. Run it and find out before the Syn2b paper repeats the claim.

---

## 7b. Closed-genome cohorts for the inversion channel

**Cost: selection is instant; the digestion + all-vs-all scales with cohort size.**

Two constraints decide which organisms are usable, and together they rule out most
of the textbook inversion systems:

- **Size.** Phase 2 measured L50 at 2,611 bp (BcgI) and 1,242 bp (four-enzyme). And
  `inverted_fraction` is length-weighted over the whole genome, so a 1 kb switch in
  a 5 Mb genome moves it by 0.02% — invisible in this channel even when its
  junctions are detected. Phase-variation invertons (Bacteroides CPS promoters
  ~200-300 bp; the Salmonella `hin`/`fljBA` flagellar switch ~1 kb, four-enzyme
  power ~0.375) are **out of scope for this channel** and belong to the junction
  channel instead.
- **Contiguity.** Per task 6, draft assemblies carry no usable orientation signal.
  Cohorts must be closed or near-closed.

Selection runs off GTDB metadata, so the accessions come out verified rather than
recalled. Check the names first — GTDB splits some NCBI species into `_A`/`_B`:

```bash
$PY $ROOT/scripts/select_closed_genomes.py \
    --metadata /lustre1/g/aos_shihuang/data/gtdb-r207/metadata/bac120_metadata_r207.tsv \
    --list-species Streptococcus
```

Then select, which also reports what is already in `genomes_all`:

```bash
$PY $ROOT/scripts/select_closed_genomes.py \
    --metadata   /lustre1/g/aos_shihuang/data/gtdb-r207/metadata/bac120_metadata_r207.tsv \
    --genome-dir /lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all \
    --species "Streptococcus pneumoniae" "Salmonella enterica" \
              "Bordetella pertussis" "Pseudomonas aeruginosa" \
    --max-contigs 3 --max-per-species 200 \
    --outdir $WORK/../closed_inversions
```

Anything not already local:

```bash
bash $ROOT/scripts/fetch_genomes_by_accession.sh \
    $WORK/../closed_inversions/accessions_to_download.txt \
    /lustre1/g/aos_shihuang/data/closed_inversions
```

Then the existing runner does the rest, on the pair file the selector wrote:

```bash
$PY $ROOT/scripts/run_syn2b_inverted_fraction.py \
    --pairs      $WORK/../closed_inversions/pairs_all_vs_all.tsv \
    --genome-dir /lustre1/g/aos_shihuang/data/closed_inversions \
    --syn2b      /lustre1/g/aos_shihuang/Syn2b/target/release/syn2b \
    --tgt-dir    $WORK/../closed_inversions/tgts \
    --out        $WORK/../closed_inversions/syn2b_inverted_fraction.tsv \
    --workers    4
```

### Target species, and why each one

GTDB names are what the selector needs. Reference strains are given only as
anchors for the locus-level work — **verify each accession before a bulk fetch**;
the selector's output is the authoritative list.

| Priority | GTDB species | Why | Event size | Reference strain |
|---|---|---|---|---|
| 1 | `Streptococcus pneumoniae` | `hsdS`/`ivr` (SpnD39III) phase-variable type I R-M system: inversions produce six alleles, changing the methylome and thus global expression and virulence | ~5-7 kb, above the BcgI floor | D39, TIGR4 |
| 1 | `Salmonella enterica` | rrn-mediated large chromosomal inversions, strongest in Typhi/Paratyphi; associated with host restriction | 10^5-10^6 bp, ideal for this channel | Typhi CT18, Typhimurium LT2 |
| 1 | `Bordetella pertussis` | >200 copies of IS481 drive large lineage-associated rearrangements; long-read structural surveys exist to check against | 10^4-10^6 bp | Tohama I |
| 2 | `Pseudomonas aeruginosa` | large inversions accumulate in cystic-fibrosis chronic infection; adaptation phenotypes (mucoidy, resistance, motility loss) are measured, and longitudinal series give a time axis | 10^4-10^6 bp | PAO1, PA14 |
| 2 | any of the above | replichore balance: inversions spanning the terminus asymmetrically reduce growth rate — a mechanism testable directly on the task-6 pairs without new phenotype data | genome-scale | — |

Filter Salmonella to Typhi/Paratyphi after selection if the full-species cohort is
too large; `ncbi_strain_identifiers` is carried through to `genomes.tsv` for that.

Two notes on scope. `Mycobacterium tuberculosis` is a poor choice despite the
phenotype data — it is structurally near-invariant. And if the 528-genome H. pylori
cagPAI cohort is Illumina drafts, the inversion channel cannot be run on it; the
cagPAI deletion analysis is unaffected, since that works from locus coordinates and
`af`, both of which survive fragmentation.

### A power-validation set that needs no phenotype

The published inverton catalogues (invertible regions called from read orientation
across human gut metagenomes) give coordinates and, more usefully, a size
distribution. Most entries fall below our floor — which is the point: the Poisson
model predicts *which* ones are detectable with no free parameters, and that
prediction can be checked against the catalogue. It converts the detection floor
from a caveat into a quantitative result. Worth doing before the phenotype cohorts,
since it needs no new genomes.

---

---

## 8. gtdb50k under FracMinHash landmarks — does the error model survive a source swap?

**Cost: one digestion pass plus one synteny pass over the same pairs. No downloads,
no alignment.** Needs Syn2b at `c10bfa3` or later on HPC (`--mode` did not exist before).

Syn2b now takes landmarks from either source (`README.md`, "Landmark sources"):

```bash
$PY $ROOT/scripts/run_syn2b_inverted_fraction.py \
    --pairs      $WORK/pairs_50k.tsv \
    --genome-dir /lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all \
    --syn2b      /lustre1/g/aos_shihuang/Syn2b/target/release/syn2b \
    --mode fracminhash --kmer 31 --scale 750 \
    --tgt-dir    $WORK/syn2b_tgts_cache_fmh750 \
    --out        $WORK/syn2b_inverted_fraction_50k_fmh750.tsv \
    --workers    16
```

**Use a separate `--tgt-dir` per mode.** The cache keys on accession alone, so
reusing the four-enzyme cache would silently feed enzyme TGTs to a FracMinHash run.
(Syn2b refuses a *mixed* comparison, but a cache hit is not a mixed comparison — it
is a run that quietly never used the mode you asked for.)

### What it answers

`scale 750` gives ~6,030 landmarks on a 4.5 Mb genome against the four-enzyme
panel's ~6,080, so this is density-matched and the two runs are directly comparable.
Three things fall out:

1. **Is the error model a property of the estimator or of the enzymes?**
   `Var(err) = 1.504*p(1-p)/m + 0.0205^2` was fitted on four-enzyme landmarks
   (R^2 0.9988 over 12 bins, out-of-sample SD(z) 1.006). If the same two constants
   come back from FracMinHash landmarks, the model is a property of the adjacency
   mathematics. If they move, they are a property of the panel — which is worth
   knowing before the paper states them as a design result.

2. **A continuous sweep of `m`, which the enzyme path cannot do.** The enzyme path
   offers four discrete densities (1, 2, 4, 16 enzymes). `--scale` is continuous, so
   the `1/m` term can be tested across a decade rather than at four points. Suggested
   rungs once the matched run is in: `--scale 250 / 750 / 2000 / 6000`.

3. **Does the sub_2 mechanism show up at scale?** Measured on one genome
   (`docs/MATH_REVIEW.md`): on a substitution ladder the four-enzyme panel leaves
   `scj_distance` 6 at 3% and 18 at 5% while FracMinHash leaves 0, because E. coli
   K-12 has 116 Hamming-1 near-duplicate landmark pairs under the panel and 0 under
   the sketch. The gtdb50k pairs are the test of whether that holds across taxa —
   compare `scj_distance` between the two runs on the same pairs, banded by ANIm.

### Expected, and what would be surprising

`breakpoints`, `inverted_fraction` and `observable_fraction` should agree closely;
verified identical on five genome-scale controls (self, 1.2 Mb rotation, 500 kb
inversion, 1/2/3/5 inversions, 40-contig shatter). A systematic gap in
`inverted_fraction` at high ANIm would be the interesting result and would need
explaining before either number is published.

Send back `syn2b_inverted_fraction_50k_fmh750.tsv`. The comparison against
`syn2b_inverted_fraction_50k.tsv` runs locally.


---

## 9. Do the two methods put breakpoints in the same *places*?

**Cost: a parse pass over files already on disk. No alignment, no digestion.**

`MATH_REVIEW.md` claims Syn2b and dnadiff place breakpoints at the same coordinates
(~1,544,640 / ~1,945,050) at a ~2790x cost difference. That rests on two hand-checked
cases. Both sides' coordinates exist for all 43k pairs already, and nothing has ever
collected them.

### Check this first — it decides whether the task is possible

`syn2b synteny` writes `synteny.junctions.tsv` next to its CSV, and
`run_syn2b_inverted_fraction.py` reads only the CSV. The pair directories are never
cleaned up, so the files *should* still be under the TGT cache:

```bash
ls $WORK/syn2b_tgts_cache/tmp_pairs | wc -l
find $WORK/syn2b_tgts_cache/tmp_pairs -name 'synteny.junctions.tsv' | wc -l
ls $WORK/out/*/dd.1coords 2>/dev/null | wc -l
```

If the first two are ~43k, this is a parse pass. If scratch was cleaned, the Syn2b
side needs re-running (cheap — the TGT cache makes it a synteny pass, not a
digestion pass); the dnadiff side cannot be recovered without re-alignment.

### The dnadiff side is *derived*, and that is the catch

`run_dnadiff_slice.sh` deletes `dd.rdiff` and `dd.qdiff` — dnadiff's own
structural-difference coordinates — to save space:

```bash
rm -f "$DD"/dd.delta ... "$DD"/dd.qdiff "$DD"/dd.rdiff "$DD"/dd.snps
```

It keeps `dd.1coords`, the 1-to-1 alignment coordinates those files are derived from,
so the boundaries can be re-derived: walk the blocks in reference order and mark
where the query stops being collinear — different query sequence (translocation),
opposite strand (inversion), or a backwards jump (relocation). Classification is by
*order*, not distance, so there is no tolerance to tune and a pure indel is correctly
not a boundary.

That is a re-derivation, not dnadiff's output, so it is checked rather than trusted:
`dd.report`'s `[Feature Estimates]` carries dnadiff's own Relocations +
Translocations + Inversions, and the script compares its derived count per pair.
**Only pairs where the two agree exactly are used** (`--all-pairs` to override). The
paper must say dd.rdiff was not retained.

```bash
$PY $ROOT/scripts/gtdb50k/collect_junction_coordinates.py \
    --pairs     $WORK/pairs_50k.tsv \
    --tgt-cache $WORK/syn2b_tgts_cache \
    --dnadiff   $WORK/out \
    --outdir    $WORK/junction_coords \
    --workers   32
```

Then, locally or on HPC:

```bash
$PY $ROOT/scripts/gtdb50k/compare_junction_positions.py \
    --coords $WORK/junction_coords/junction_coordinates.tsv \
    --truth  $WORK/high_ani_truth.tsv \
    --syn2b  $WORK/syn2b_inverted_fraction_50k.tsv \
    --out    $WORK/junction_coords/position_agreement.tsv
```

### Reading the result

The two sets are matched **one-to-one**, greedily by increasing distance — nearest
neighbour in both directions would let one Syn2b junction satisfy three dnadiff
boundaries and read as three successes. Greedy is not the optimal assignment, but it
can only inflate matched distances, so the agreement reported is a lower bound.

**The expected answer is set by landmark spacing, not by the algorithm.** Syn2b
reports the left landmark of the broken adjacency, so its error is bounded by the gap
to the next landmark. Measured on a closed E. coli K-12 control with three known
200 kb inversions:

| landmarks | spacing | median error | max error |
|---|---|---|---|
| BcgI, 2,872 | 1,582 bp | 614 bp | 3,671 bp |
| 4-enzyme, 6,079 | 747 bp | 273 bp | 3,671 bp |
| FracMinHash s=750, 6,034 | 753 bp | 456 bp | 1,031 bp |
| FracMinHash s=200, 22,708 | 200 bp | 248 bp | 783 bp |
| FracMinHash s=50, 90,394 | 50 bp | 44 bp | 71 bp |

So a median matched distance near the panel's spacing is the *expected* result and
the correct claim. A distance much larger than spacing is the finding that would need
explaining. Note the enzyme path's max error does not improve from BcgI to the
four-enzyme panel — restriction sites cluster on their motifs, so adding enzymes does
not fill a gap that has none — while FracMinHash at matched density already reaches
1,031 bp because its spacing is Poisson-uniform.

Send back `junction_coordinates.tsv` and `position_agreement.tsv`.


## Not needed any more

**Resampling pairs at >=97% ANIm.** This was on the list because held_out_50k has
n = 2 there. The high-ANI set already covers it: 3,826 pairs at >=97% measured ANIm,
with slope 1.0063, r = 0.9960 and SD 0.0135. Combined with held_out_50k's 80–97%
coverage, the ANIm axis is continuous from 80% to 100% with no gap. See the banded
tables in `inverted_fraction_comparison_report.md`.
