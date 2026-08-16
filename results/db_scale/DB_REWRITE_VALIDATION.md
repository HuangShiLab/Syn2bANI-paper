# Database-path rewrite validation: screen + chain-restricted MLE for `dist`/`search`/`triangle`/`db`

Date: 2026-08-17. Author: kimi-code agent (agent-12). No git commits made — diff under review.
Companion to [DB_SCALE_BENCHMARK.md](DB_SCALE_BENCHMARK.md), which documented the
legacy-path failures this rewrite fixes.

## What changed (syn2bani working tree on top of 15da386)

All database-scale subcommands now run a two-stage pipeline on the validated
estimator, sharing one implementation in `src/cli/compare.rs`:

- **Stage 1 — screen** (`src/core/screen.rs`): per enzyme, each tag is reduced
  to one strand-canonical key (its centred 18 bp packed window); keys are
  sorted/deduplicated per genome, and a pair passes iff shared keys ≥ 3 AND
  shared/min(smaller key set) ≥ 0.001. Both conditions required (see
  calibration below). Tunable: `--screen-window`, `--screen-min-shared`,
  `--screen-min-containment`.
- **Stage 2 — refine**: survivors go through `chain_ani::compute` via the
  identical call and row formatter `ani` uses (`compare::refine_pair`,
  `compare::ani_row`), so numbers are byte-identical to an `ani` run on the
  same pairs. Output = `ani` TSV columns plus a trailing `flag`
  (`ok`/`INCONSISTENT`/`BELOW_DETECTION`).
- Optional second-tier gate `--refine-min-approx` (crude containment-ANI;
  default 0 = off) exists to bound refine calls; it was **not needed** —
  measured refine cost is acceptable (below).

Per subcommand:

- `dist`: `--ql/--rl` added (positional `<queries...> <reference>` kept,
  ani-style). Every genome is digested/loaded **once** (the legacy path
  re-digested every reference per query). `--enzymes` default now
  `BcgI,AlfI,AloI,FalI`. The legacy flags `--min-af`, `--mash-ani`,
  `--raw-features`, `--structural`, `--multi-enzyme`, `-e` are **removed**
  (breaking change; the GBRT feature-dump mode that training scripts under
  `scripts/run_*` used is gone).
- `search`: queries (FASTA or `.s2ba`, positional or `--ql`) against a sketch
  directory (positional) or `--rl` list of sketch paths. No more hardcoded
  BcgI; the DB sketches' recorded panel is honoured (v1 sketches without an
  enzyme table are refused with an explicit error). Hits filtered by
  `--min-ani` (default 0.8) on `ani_gated`, reported best-first per query.
- `triangle`: positional list or `--ql`; `--edge-list` emits ani-style rows
  for refined pairs with a finite gated estimate; matrix mode writes `NaN`
  (not `0.0000`) for screened-out/below-floor pairs. No more hardcoded BcgI.
- `sketch` / `db build`: default panel changed from BcgI-only to
  `BcgI,AlfI,AloI,FalI` (**breaking change**, documented in README; `.s2ba`
  files record their enzyme table and old sketches stay readable).
- `db add`: digests with the panel recorded in the existing DB (was hardcoded
  BcgI); empty DB falls back to the default panel.
- `db search`: same screen→refine pipeline for sketch-dir × sketch-dir.
- The v7 `TagMatcher`/`AniCalculator`/GBRT path is no longer used by any of
  these subcommands. The modules remain in the crate (`core/gbrt.rs` etc.)
  but nothing on the db path calls them; `ani --calibrate` uses the separate
  linear model in `core/calibration.rs` and is unaffected (verified:
  `ani --calibrate` adds the `ani_cal` column as before).

`src/cli/ani.rs` was refactored, not rewritten: `load_sketch` is now
`pub(crate)`, the TSV header/row construction moved verbatim into
`compare::ani_header`/`compare::ani_row`, and `ani` calls them — so `dist`
and `ani` cannot drift apart.

## Unit/local validation (macOS, cargo test --release)

- 91 tests pass (80 lib incl. 7 new `core::screen` tests: strand invariance,
  window sensitivity, per-enzyme isolation, threshold logic, approx
  monotonicity; 11 integration).
- realbench (14 complete enteric chromosomes, `prototype/realbench/genomes`):
  `dist --ql --rl` all-vs-all is **byte-identical to `ani`** on all 196 pairs
  (max |Δani_gated| = 0.0); self-comparisons 99.9999; 196 pairs in 0.7 s
  at 8 threads. `search` on `.s2ba` sketches matches `ani` on FASTA exactly
  (max |Δ| = 0 over 56 hits). `db build`/`db add`/`db search`/`triangle`
  (edge-list and matrix) smoke-tested; `db add` inherits the DB panel.

## Screen calibration (HPC, GTDB-R207)

Screen statistics were dumped for window widths 13/15/17/18 on two sets:
(a) the cross product of the 996 genomes in `lists/accuracy_pairs.tsv` —
500 validated true pairs (ANI 80–100, v8 MLE truth) for recall;
(b) the n=500 random GTDB subset cross product (250k pairs, skani finds zero
≥80% pairs) for selectivity. Raw dumps on HPC at
`results/db_scale/rewrite/screen_{acc,n500}_w{13,15,16,17,18}.tsv`.

Why 18 bp: at 11–13 bp generic bacterial background homology (shared genes at
~65–75% identity) makes key containment useless as a screen (random-pair
median containment 0.47 at W=11, 0.15 at W=13). At 18 bp the bulk of random
pairs shares 0–2 keys.

| gate (W=18) | FRR on 500 true pairs | random-GTDB rejection |
|---|---|---|
| shared≥3 AND cont≥0.001 | **0/500** | 82.7% |
| shared≥4 AND cont≥0.0015 | 0/500 | 88.4% |
| shared≥4 AND cont≥0.002 | 1/500 (0.2%) | 91.2% |

Shipped default: **W=18, shared≥3 AND cont≥0.001** (recall-first). Margins:
in the FRR-critical 80–85% band the weakest true pair has shared=29,
cont=0.0061 — 10×/6× above the floors. The pairs closest to the floor are
95–100%-band calls with kilobase-scale overlap (see caveats).

The legacy exact-tag `min_af=0.1` screen rejected **94.4%** of true ≥80%
pairs; the new screen rejects **0%** of them.

## HPC scale check (partition amd, 32 threads, /usr/bin/time -v)

Code deployed to the clean clone `/lustre1/g/aos_shihuang/Syn2bANI-bench` via
tar-over-ssh (the dirty main checkout was not touched). The 5000-genome DB was
re-sketched with the new default panel (`s2ba_n5000_new`, 479 MB store, 30 s,
1.2 GB RSS — consistent with the old explicit-`--enzymes` build at 491 MB).

| phase | rewritten | legacy v7 | skani 0.3.2 |
|---|---|---|---|
| triangle n=500 | **8.5 s / 0.88 GB** | 107.0/99.9 s / 2.41 GB | 1.0–1.6 s / 0.62 GB |
| triangle n=2000 | **86.6 s / 1.84 GB** | 893.3 s / 6.80 GB | 1.71 s / 1.83 GB |
| triangle n=5000 | **406.7 s / 4.11 GB** | not attempted (~2.5 h projected, meaningless output) | not measured |
| search 100×5000 | **22.0 s / 3.2 GB** | 558.3/550.9 s / 1.57 GB | 3.14/3.16 s / 1.15 GB |

Screen pass rates: 21,595/124,750 (17.3%) at n=500; 305,996/1,999,000 (15.3%)
at n=2000; 1,570,949/12,497,500 (12.6%) at n=5000 — consistent with the 83%
calibrated rejection. Survivors are overwhelmingly same-family background
pairs that the estimator then correctly measures or marks `BELOW_DETECTION`.

### triangle n=2000 output reconciliation

Against a run with the screen forced to pass everything (1,999,000 pairs
refined): the final screen retains **843/844** pairs that the estimator
reports a finite ANI for, with **max |Δani_gated| = 0.0** on retained pairs
(same code path). The one lost pair (JADJDU010000001.1 × DRLG01000001.1,
gated 98.44) has af_query = af_reference = **0.0000** — a tiny shared-island
call, INCONSISTENT-flagged. Measured false-reject rate on
estimator-reportable pairs: **1/844 = 0.12%** (target ≤2%).

Of the 843 reported rows: min gated ANI 83.26, median 93.63; 758/843 have
max(AF) < 0.15, i.e. skani's default AF≥15 rule would not report them at all
(skani's n=2000 matrix is entirely NA off-diagonal). The legacy triangle
reported 2,107 hits of which **0** were confirmed; the rewrite reports only
estimator-backed rows with honest AF and flag columns.

### search 100 queries × 5000 DB

184 hits at `--min-ani 0.8`. All **12/12** pairs skani reports are found, and
on the 7 intersection pairs with independent MLE truth the new search matches
that truth **exactly** (91.09/89.34/94.47/93.73/93.09/94.57/95.11). The
legacy search underestimated these same pairs by 8–11 points and returned
0.90 for a self-hit; the rewrite returns 99.9999 for self-hits. 13/184 hits
have max(AF) ≥ 0.15 (skani-comparable); the remaining 171 are low-AF calls
reported with their true AF.

### triangle n=5000

12,497,500 pairs screened in 406.7 s wall (4.1 GB RSS): 12.6% passed the
screen, and the edge list carries **3,410 estimator-backed rows** (banding of
the gated estimate: 32 rows 80–85, 735 at 85–90, 1,392 at 90–95, 1,251 ≥95;
flags: 1,719 ok, 1,382 INCONSISTENT, 309 BELOW_DETECTION). No NaN flood, no
spurious 0.0000 rows — contrast with the legacy triangle, whose every ≥80%
call was spurious.

### Bugs caught during this validation

- The first screen (11 bp window, OR-gate) passed 100% of random GTDB pairs —
  background homology makes short-window containment useless; fixed by the
  calibrated W=18 AND-gate above.
- The initial `--min-ani` implementation compared a fraction flag against a
  percent-scale value, so `-m 0.95` filtered at 0.95%. Caught by a local
  smoke test, fixed (`ani_gated >= min_ani`, fraction vs fraction) in
  `dist`/`search`/`db search`, and the HPC search was rerun: the hit set is
  unchanged (all refined hits were ≥82.9 anyway; v3 and v4 outputs are
  byte-identical).

## Caveats / honest limitations

- **Runtime vs skani**: ~50× slower on triangle n=2000 (87 s vs 1.7 s) and
  ~7× on search (22 s vs 3.1 s); triangle n=5000 completes in under 7
  minutes. The refine stage is the honest cost of running a chaining MLE per
  surviving pair; skani's regressed k-mer estimate is cheaper.
  Seconds-to-minutes scale holds through n=5000.
- **Low-AF island pairs**: the estimator reports ANI over chained regions;
  pairs sharing only a small (e.g. mobile-element) island get high ANI at
  AF≈0.001–0.01 with `ok` or `INCONSISTENT` flags. They dominate the
  triangle/search hit lists on dereplicated GTDB. skani hides this class via
  AF≥15. Consumers should filter on `af_query`/`af_reference`/`flag`; the
  screen is deliberately loose enough to keep them (one such pair was still
  lost at n=2000 — AF exactly 0).
- **`--refine-min-approx`** exists as an emergency brake for huge,
  relatedness-rich inputs but is off by default and unneeded at n≤5000.
- Screen thresholds are calibrated on the 4-enzyme panel; a different panel
  changes key counts and should be recalibrated (flags are exposed).

## Files

- Code (syn2bani repo, uncommitted): `src/core/screen.rs` (new),
  `src/cli/compare.rs` (new), `src/cli/{ani,dist,search,triangle,db,sketch,mod}.rs`,
  `src/core/mod.rs`, `src/main.rs`, `README.md`.
- This repo: `results/db_scale/rewrite/{rewrite_scaling*.tsv,
  triangle_n500_v3.tsv, triangle_n2000_v3.tsv, triangle_n5000_v3.tsv,
  search_100x5000_v3.tsv, search_100x5000_v4.tsv}`.
- HPC only (lustre, `results/db_scale/rewrite/`): screen calibration dumps
  `screen_{acc,n500}_w*.tsv`, pass-all-screen triangle outputs
  `triangle_n{500,2000}.tsv`, `search_100x5000.tsv` (v1), logs, SLURM
  scripts; sketch store `sketches/s2ba_n5000_new`.
