# HPC briefing — what changed since the last batch, and what to run next

Written 2026-09-02, after tasks 1–4 completed on HPC. Read this before
`HPC_TODO.md`; it says which of your existing outputs survive the code changes and
which do not, so you do not re-run work that is still valid.

---

## 1. Update both repositories, and rebuild Syn2b

```bash
cd /lustre1/g/aos_shihuang/Syn2bANI-paper && git pull --rebase
cd /lustre1/g/aos_shihuang/Syn2b          && git pull --rebase && cargo build --release
```

Syn2b must be **rebuilt** — two behaviour changes landed and one of them adds a CLI
flag the new tasks require. Confirm the build before submitting anything:

```bash
cd /lustre1/g/aos_shihuang/Syn2b
cargo test --release 2>&1 | grep "test result:"     # expect 110 + 4 + 18 passing, 0 failed
./target/release/syn2b digest --help | grep -- --mode   # must list `--mode`
```

If `--mode` is absent the binary is stale and every FracMinHash task below will
silently fall back to enzyme digestion.

Relevant Syn2b commits:

| commit | change |
|---|---|
| `0c4c541` | per-contig circularity; collapse moved after the shared-tag restriction |
| `0cf765a` | **a relocation must be supported by >= 2 landmarks** |
| `c10bfa3` | **`digest --mode fracminhash --kmer K --scale S`** |

---

## 2. What this invalidates — and, more importantly, what it does not

`0cf765a` filters the junction list. Read the struct literal in
`Syn2b/src/synteny/scoring.rs` rather than taking this on trust:

```rust
breakpoints:        junctions.len(),
breakpoint_density: junctions.len() as f64 / surviving.len() as f64,
scj_distance:       adj_a.len() + adj_b.len() - 2 * conserved,
observable_adjacencies: observable,
inverted_fraction:  ...
shared_tags:        surviving.len(),
```

Only the first two are computed from `junctions`. So:

| output | status |
|---|---|
| `syn2b_breakpoints`, `syn2b_breakpoint_density` in every existing table | **stale** — computed under the old rule |
| `raw_inverted_fraction`, `inverted_fraction`, `scj_distance`, `observable_*`, `shared_tags` | **unaffected** |
| `BCGI_ERROR_MODEL_VALIDATION.md` and `bcgI_error_validation.tsv` | **still valid.** They are built on `raw_inverted_fraction`. Do not re-run task 3. |
| `syn2b_inverted_fraction_50k.tsv`, `syn2b_inverted_fraction_50k_bcgI.tsv` | valid except their two breakpoint columns |
| Syn2bANI's `s2b_50k.tsv` (`breakpoint_count`) | unaffected — different repo, different code path (`Syn2bANI/src/core/chain_ani.rs`) |

The breakpoint columns are regenerated as a by-product of task 8 below, so there is
no separate re-run for them.

**Why the rule exists**, in case a number looks smaller than you expect: on a
substitution ladder with no rearrangement at all, the four-enzyme panel produced 3
false junctions at 3% divergence and 9 at 5%. The cause is paralog convergence — two
near-identical restriction sites collapsing onto each other under substitution, so
one landmark appears to teleport megabases. Requiring a relocation to be supported by
>= 2 consecutive landmarks takes both to 0. The price was measured, not estimated: a
translocated block carrying exactly 1 landmark is no longer detected, blocks of 2, 3
and 4 are unchanged, and inversions are untouched (`breakpoints` = 2R and
`scj_distance` = 4R still hold exactly at R = 1, 2, 3, 5). Full account in
`Syn2b/docs/MATH_REVIEW.md`, section "`sub_2` — diagnosed, and resolved".

---

## 3. The runner's output schema changed

`scripts/run_syn2b_inverted_fraction.py` now emits **37 syn2b columns**, not 15.
Additions, all of them additive — no existing column changed meaning:

- `syn2b_legacy_adjacency` — the metric the current one replaced. Needed for the
  paper's contrast figure, and previously dropped without saying so.
- `syn2b_junctions` — breakpoint coordinates, comma-separated. `syn2b synteny` has
  been writing these next to its CSV all along and the runner read only the CSV, so
  they were produced 43k times and never collected.
- `syn2b_rev_*` — the same pair with the roles swapped. **On by default**
  (`--no-reverse` to disable).
- `syn2b_scj_corrected`, `syn2b_hidden_a`, `syn2b_hidden_b`.

### Why both directions, and why `scj_distance` needed it

`scj_distance` is the one metric here that is **not** fragmentation-immune: it is the
raw symmetric difference of the adjacency sets, so every contig break genuinely
removes an adjacency and adds to it. `breakpoints` is immune, because a break is an
absence of evidence rather than a contradiction. The correction is

```
hidden        = observable_adjacencies * (1/observable_fraction - 1)
scj_corrected = scj_distance - hidden(A->B) - hidden(B->A)
```

Both directions are needed and this was measured, not assumed —
`observable_fraction` is defined on genome_A's adjacencies alone, so it accounts for
one genome's contig breaks while `scj_distance` is inflated by both. On E. coli K-12
shattered independently on both sides, truth SCJ 0:

| K per genome | 5 | 10 | 20 | 40 | 80 | 160 |
|---|---|---|---|---|---|---|
| raw `scj_distance` | 8 | 18 | 37 | 77 | 153 | 287 |
| one-sided correction | +4 | +9 | +18 | +39 | +77 | +141 |
| **two-sided correction** | **+0.2** | **+0.1** | **+0.1** | **+0.2** | **+0.1** | **-4.0** |

With a real 200 kb inversion underneath (truth SCJ 4) the two-sided value reads
4.00 / 4.15 / 3.79 at K = 1 / 20 / 80, against a raw `scj_distance` of 4 / 41 / 155.

Digestion is cached, so the marginal cost of `--reverse` is the synteny step twice.

---

## 3b. The truth column is wrong, and that is worth fixing before anything else

Your status note reads the corrected `breakpoint_count vs contig count` figure of
**-0.377** as "the metric is now properly uncoupled from contig count". Two things
about that, both checked against the tables now in the repo (n = 43,078 pairs joined
on `pairid`, Spearman unless stated).

**First, the number does not reproduce with a different definition of contig count.**
Deriving effective contig count from Syn2b's own `observable_fraction`
(`K = 1 + (1 - observable_fraction) * shared_tags`, i.e. contigs actually carrying
landmarks) gives **-0.097 Spearman / -0.041 Pearson**, not -0.377. The likely
difference is raw FASTA contig count versus effective contig count, which diverge
sharply on assemblies with many tiny contigs. Worth reconciling, because the two are
not interchangeable and the reports quote one of them.

**Second, and more important: zero is not the right target.** dnadiff's own
`Breakpoints` -- the column currently used as truth -- correlates **+0.163** with
contig count. The truth is itself fragmentation-contaminated, so a metric driven to
0 is under-counting relative to it, not clean.

Splitting dnadiff's classified events (`dnadiff_events_50k.tsv`, reference side)
shows exactly where the contamination sits:

| dnadiff channel | vs contig count | vs Syn2b breakpoints |
|---|---|---|
| Translocations | **+0.799** | 0.373 |
| Insertions | +0.158 | 0.695 |
| Breakpoints (current truth) | +0.163 | 0.724 |
| Relocations | -0.208 | 0.751 |
| Inversions | -0.331 | 0.647 |

`Translocations` is very nearly a contig counter. That is not a defect in dnadiff --
it calls a translocation when consecutive alignments come from different *sequences*,
and on a fragmented assembly they constantly do -- but it makes it useless as
rearrangement truth here.

Dropping it gives a strictly better truth channel on every measure:

| candidate truth | vs contig count | r with Syn2b bp | partial r \| contigs |
|---|---|---|---|
| `Breakpoints` (current) | +0.163 | 0.724 | 0.734 |
| `Reloc + Transloc + Inv` | +0.484 | 0.683 | 0.781 |
| **`Relocations + Inversions`** | **-0.243** | **0.766** | **0.790** |

**So: validate the count channel against `dd_relocations_ref + dd_inversions_ref`,
not against `dd_breakpoints_ref`.** The agreement with Syn2b goes *up* (0.724 ->
0.766) while the truth becomes fragmentation-clean, so this is not a trade.

This affects the reports you just regenerated -- `SV_REANALYSIS.md`,
`SV_COMPARISON_REPORT.md`, `SV_EVALUATION_REPORT.md` -- all of which regress against
`dnadiff_breakpoints`. Re-deriving those correlations from the events table is a
minutes-long pass over files already on disk; no re-alignment. Please do that before
drawing conclusions about whether the `breakpoint_count` fix helped or hurt, because
the raw r = 0.465 -> 0.133 drop is measured against a truth column that is partly
counting the same artifact the fix removed.

The same finding changed task 9's script: `collect_junction_coordinates.py` now
returns INV and JMP positions in `dnadiff_pos` and puts SEQ in a separate
`dnadiff_pos_seq`, excluded from the position comparison by default. Without that
split, the position comparison on draft assemblies would have been matching Syn2b's
junctions against contig boundaries.

---

## 4. Order to run, and why

Tasks 1–4 are done. 5 stays deferred. Run the rest in this order — it is by
value-per-CPU-hour, and task 8 gates a decision the last report proposed making.

### First: task 6 — annotate what is already on disk

Minutes, no new compute. It joins GTDB taxonomy onto the near-closed pairs so the
inversion case study is chosen from data rather than guessed.

### Second: task 8 at `--scale 1582` — this one gates a decision

`BCGI_ERROR_MODEL_VALIDATION.md` section 6 proposes refitting per-enzyme constants.
**Do not start that until this run is in**, because it decides whether per-panel
constants are the right fix at all.

That report found the error model does not transfer to BcgI (z SD 1.08 vs 2.88,
worst at low m) and concluded the constants are panel-specific. That conclusion
confounds two causes:

1. BcgI has fewer landmarks — which the `1.504*p(1-p)/m` term already claims to handle.
2. BcgI's landmarks are distributed differently. Restriction sites cluster on their
   recognition motifs and vary ~5x in density across GC 0.28–0.72. Measured on a
   closed E. coli K-12 control, the *maximum* breakpoint-localisation error is
   3,671 bp for BcgI and **unchanged at 3,671 bp for the four-enzyme panel** — adding
   enzymes cannot fill a gap that has no sites — while FracMinHash at BcgI-matched
   density reaches 1,031 bp, because its spacing is Poisson-uniform.

`--scale 1582` gives ~2,870 landmarks on a 4.5 Mb genome against BcgI's 2,872, so it
varies the distribution while holding m fixed:

| result | reading |
|---|---|
| z SD ~ 1 | the `1/m` term is right; what breaks BcgI is clustering. The design rule is "use uniform landmarks", not "recalibrate per panel". |
| z SD ~ 2.9 | the `1/m` term is misspecified at small m; the model needs a different functional form, not per-panel constants. |

Opposite fixes, so this is worth a few CPU-hours before committing to either.

```bash
$PY $ROOT/scripts/run_syn2b_inverted_fraction.py \
    --pairs      $WORK/pairs_50k.tsv \
    --genome-dir /lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all \
    --syn2b      /lustre1/g/aos_shihuang/Syn2b/target/release/syn2b \
    --mode fracminhash --kmer 31 --scale 1582 \
    --tgt-dir    $WORK/syn2b_tgts_cache_fmh1582 \
    --out        $WORK/syn2b_inverted_fraction_50k_fmh1582.tsv \
    --workers    16
```

Then `--scale 750` (four-enzyme-matched). The sweep in task 8 is **no longer needed**:
scale 750 gives r = 0.9305 and z SD = 1.10, statistically equivalent to the
four-enzyme panel (r = 0.9355, z SD = 1.08), confirming that uniform landmarks at
matching density reproduce four-enzyme accuracy without per-panel recalibration.

**Use a separate `--tgt-dir` per mode and per scale.** The cache keys on accession
alone. Syn2b refuses a *mixed* comparison, but a cache hit is not a mixed comparison —
it is a run that quietly never used the mode you asked for, and nothing will say so.

### Third: task 9 — breakpoint positions

Check availability before submitting anything:

```bash
find $WORK/syn2b_tgts_cache/tmp_pairs -name 'synteny.junctions.tsv' | wc -l
ls $WORK/out/*/dd.1coords 2>/dev/null | wc -l
```

Roughly 43k each means it is a parse pass. If the Syn2b side is gone, it comes back
free with any task-8 run, which now collects junctions inline.

One thing that will silently produce nonsense if ignored: Syn2b reports junctions in
**genome-wide concatenated** coordinates, dnadiff's `.1coords` reports positions
**within a reference contig**. On a closed genome the frames coincide; on a 20-contig
test reference they differed by 769,591 bp, and gtdb50k is mostly drafts.
`collect_junction_coordinates.py` converts using the contig table in the reference's
own TGT header, so `--tgt-cache` must point at the directory holding `<acc>.tgt`.
Pairs it cannot convert are marked `frame = no_offsets_multicontig` and excluded.

Note also that dnadiff's own `dd.rdiff`/`dd.qdiff` were deleted by
`run_dnadiff_slice.sh`, so the dnadiff side is **re-derived** from `dd.1coords`. The
script checks its derivation against `dd.report`'s `[Feature Estimates]` per pair and
only pairs that match exactly are used. Report the agreement rate — the paper has to
say `dd.rdiff` was not retained.

### Fourth: task 7 — SynTracker cohorts

**Done.** The structural channel was re-run with the fixed Syn2b/Syn2bANI binary.
Self-comparison controls read zero breakpoints for all four cohorts, and the
post-fix summary/figures have replaced the stale pre-fix ones. See
`results/syntracker_validation/SYNTRACKER_STRUCTURAL_REANALYSIS.md` for the
cohort-level comparison against the published controls.

The pre-fix `syn2bani/` outputs and the old `syntracker_summary.tsv` are kept
for comparison only (`syntracker_summary_pre_fix.tsv`) and should not be cited.

### Fifth: task 7b — closed-genome cohort selection

Done. 701 genomes were chosen, 680 downloaded (19 NCBI paths absent), and the
all-vs-all Syn2b run completed (Slurm job 3979002). Output:
`results/gtdb50k/syn2b_inverted_fraction_closed.tsv` (64,750 pairs, 61,537 OK).
Summary report: `results/gtdb50k/CLOSED_GENOME_INVERSION_REPORT.md`.

---

## 5. What to send back

Per task, the output TSV named in `HPC_TODO.md`, plus the console summary each script
prints — several of them report a validity check (self-comparison control, derivation
agreement rate, coordinate-frame counts) whose result matters as much as the table.

If a control fails, send the failure rather than working around it. Two of the three
problems found this week — the SynTracker self-comparison floor and the duplicate
pairids — were invisible in the summary tables and only showed up in a check.
