# Database-scale benchmark: syn2bani `sketch`/`dist`/`search`/`triangle` vs skani

Date: 2026-08-14. Author: automated benchmark run (kimi-code agent), no git commits made.

## Verdict (one paragraph)

`syn2bani sketch` is healthy: linear scaling, ~3.5x smaller sketch store than skani,
~8x slower wall time but seconds-to-minutes at n=5000. **The comparison
subcommands (`dist`, `search`, `triangle`) are not production-ready**: none of them
uses the validated chain-restricted MLE estimator of `ani` — they all call the
legacy v7 `TagMatcher`/`AniCalculator` (GBRT-debias) path that `struct.rs` already
flags as untrustworthy. Consequences measured below: `triangle` reports **only
spurious hits** on a GTDB sample where skani finds zero >=80% ANI pairs (0/316 and
0/2107 confirmed), `search` underestimates genuine hits by ~8-11 ANI points and
floods the output with threshold-edge false positives, and `dist` with default
settings screens out 94.4% of true >=80% ANI pairs (though it is accurate to
~0.9 ANI points on the few pairs it does report).

## Environment

- HPC: hpc2021.hku.hk. Phases n<=500 (sketch, triangle) ran on the **login node**
  (64 cores, shared — timings noisy, one 36 s I/O outlier at n=100); n>=2000
  phases ran as SLURM job 3907622 on partition `amd`, node GPA-2-4
  (AMD EPYC 7742, 128 cores allocated 32, 180 GB).
- Threads = 32 everywhere (`-p -t 32` / `skani -t 32`).
- syn2bani 0.1.0, commit **98177dc** ("ani: gated estimator and recalibrated
  consistency flag"), `cargo build --release`, rustc 1.97.1. Built in a clean
  clone at `/lustre1/g/aos_shihuang/Syn2bANI-bench` because the main checkout
  `/lustre1/g/aos_shihuang/Syn2bANI` has uncommitted changes to
  `src/cli/ani.rs`, `src/core/chain_ani.rs`, `src/core/calibration.rs`,
  `ALGORITHM_MLE.md` that blocked `git pull --ff-only`. The dirty tree was left
  untouched.
- skani 0.3.2 (conda env `syn2bani`, `/home/shihuang/.conda/envs/syn2bani/bin/skani`).
  skani ran in its default "learned ANI mode" (regression-adjusted ANI).
- Genome pool: 65,703 GTDB-R207 assemblies,
  `/lustre1/g/aos_shihuang/data/gtdb-r207/genomes_all/*.fna`.
- Subsets: nested random samples (fixed seed, `shuf --random-source=<(yes 42)`):
  n100 < n500 < n2000 < n5000; 100 held-out queries = rows 5001-5100 of the same
  permutation (disjoint from n5000). Lists in `lists/`.
- Timing: `/usr/bin/time -v` (wall, max RSS). Reps: 3 for n=100/500 sketch,
  2 for n=500 triangle / n=2000 sketch / search, 1 for n=5000 sketch and
  n=2000 triangle.
- syn2bani sketch panel: `--enzymes BcgI,AlfI,AloI,FalI` (the validated `ani`
  default panel) for all sketch/search work.

Exact commands: see `scripts/` (`bench_sketch.sh`, `bench_triangle.sh`,
`bench_search.sh`, `bench_accuracy.sh`, `slurm_heavy.sh`, `common.sh`).

## Subcommand semantics (from source at 98177dc + smoke tests)

| subcommand | estimator | enzyme panel | notes |
|---|---|---|---|
| `ani` | chain-restricted MLE (`chain_ani::compute`) — validated | `-e BcgI,AlfI,AloI,FalI` default | reference path |
| `sketch` | n/a (digest only) | default `-e BcgI` only; `--enzymes` accepts full panel | records enzyme table in `.s2ba` |
| `dist` | v7 `TagMatcher` + `AniCalculator` | default `-e AloI,BslFI`; `--enzymes` accepts panel; multi-enzyme disables near-match | two greedy positionals: **all-but-last arg = queries, last arg = single reference**; no `--ql/--rl`; re-parses+re-digests every reference per query (O(QxR) I/O); default output = chained-kmer mash-like ANI, `--mash-ani` = GBRT v7; coarse screen `min_af=0.1` default; `mash_calibration_offset` exists in config, always 0.0 |
| `search` | v7 path, GBRT v7 output | **query hardcoded to BcgI**; DB tags lumped with `enzyme="unknown"`, `contig_id=0` (DB panel recorded but ignored) | no enzyme flag exists; `min_ani=0.8` default; loads entire DB into RAM |
| `triangle` | v7 path, GBRT v7 output | **hardcoded BcgI, no enzyme flag** | all-pairs in memory; writes `0.0000` (not NaN) for unrelated pairs; keeps full sequences in RAM |
| `db build` | digest | single `-e` or `--multi-enzyme` (all registry enzymes) — **cannot take the 4-enzyme panel** | delegates to `sketch` |
| `db add` | digest | **hardcoded BcgI** — silently inconsistent with a panel-built DB | |
| `db search` | v7 path, GBRT v7 | panel-agnostic (lumps tags) | sketch-vs-sketch version of `search` |

**Cannot take the 4-enzyme panel:** `search`, `triangle`, `db build` (panel of 4),
`db add`. Only `sketch` and `dist` accept `--enzymes BcgI,AlfI,AloI,FalI`.

## 1. Sketch scaling (syn2bani `--enzymes BcgI,AlfI,AloI,FalI` vs `skani sketch`)

| n | tool | wall s (reps) | peak RSS GB | store MB |
|---|---|---|---|---|
| 100 | syn2bani | 9.4 / 36.2 / 5.2 (login) | 0.28 | 11.8 |
| 100 | skani | 2.6 / 8.8 / 6.1 (login) | 0.20 | 40.0 |
| 500 | syn2bani | 9.9 / 2.2 / 3.5 (login) | 0.52 | 52.3 |
| 500 | skani | 11.6 / 3.8 / 5.0 (login) | 0.53 | 181.2 |
| 2000 | syn2bani | 13.0 / 3.3 | 0.77 | 201.4 |
| 2000 | skani | 1.0 / 1.0 | 0.52 | 705.5 |
| 5000 | syn2bani | 20.3 (single) | 1.22 | 491.1 |
| 5000 | skani | 2.5 (single) | 1.10 | 1745.3 |

syn2bani sketch scales linearly (~0.1 MB/genome store, 3.5x smaller than skani),
~8x slower than skani but absolutely fast. Login-node reps are I/O-noisy
(lustre); the compute-node numbers are cleaner.

## 2. All-vs-all (`syn2bani triangle --edge-list` vs `skani triangle`)

| n | tool | wall s | peak RSS GB | output |
|---|---|---|---|---|
| 100 | syn2bani | 11.4 | 1.10 | 4,950 edges |
| 100 | skani | 0.43 | 0.18 | 100x100 matrix |
| 500 | syn2bani | 107.0 / 99.9 | 2.41 | 124,750 edges |
| 500 | skani | 1.6 / 1.0 | 0.62 | 500x500 matrix |
| 2000 | syn2bani | 893.3 (single) | 6.80 | 1,999,000 edges |
| 2000 | skani | 1.71 (single) | 1.83 | 2000x2000 matrix |

syn2bani triangle scales as expected for all-pairs (O(n^2) pair evaluations with
the v7 matcher, ~0.5 ms/pair at 32 threads) and is ~100-500x slower than skani.
RSS grows ~linearly (all digests + full sequences held in RAM). No crash at
n=2000; n=5000 not attempted (~2.5 h projected, and output is meaningless — see
below).

**Quality: every reported triangle hit is spurious.** skani finds **zero** pairs
>=80% ANI in both the n=500 and n=2000 samples (expected — GTDB-R207 genomes_all
is species-dereplicated). syn2bani triangle reports ani>0 for 316/124,750 pairs
at n=500 and 2,107/1,999,000 at n=2000, of which **0 are confirmed by skani**.
Distribution at n=2000: 1,714 pairs at 0.80-0.83, 249 at 0.83-0.87, 101 at
0.87-0.90, 32 at 0.90-0.95, **11 at >=0.95** — median evidence is 14 shared tags
(minimum-reporting threshold is 10). The v7 matcher converts a handful of
spurious shared BcgI tags into ANI 0.80-0.96.

## 3. Search (100 held-out queries x 5000-genome sketch DB)

| tool | wall s (2 reps) | peak RSS GB | reported pairs |
|---|---|---|---|
| syn2bani `search --min-ani 0.8` | 558.3 / 550.9 | 1.57 | 150 (34 queries) |
| skani `search` (default min-af 15) | 3.14 / 3.16 | 1.15 | 12 (10 queries) |

Note the recall gates differ (syn2bani filters on ANI>=80; skani on AF>=15%).
syn2bani search digests each query **BcgI-only** and matches against the 4-enzyme
DB with tags lumped enzyme-agnostic — so it is effectively a BcgI-only search
against a 4x larger tag set, with the v7 matcher. It also rebuilds every DB
tag-set per query (O(Q x DB) unpack work), hence ~9 min vs skani's 3 s.

Quality (intersection of reported pairs, `out/search_intersection_ani.tsv` and
`analyze_search.py` output):

- 10 (query,ref) pairs reported by both tools: Pearson r = 0.31, mean |dANI| =
  6.0 points; top-hit genome agrees for 6/8 queries with any hit in both.
- Independent `ani` (validated MLE) re-run on 7 of those intersection pairs:
  truth 89.3-95.1; skani 87.0-92.9 (close); **syn2bani search 80.6-84.1 —
  systematically 8-11 points low**, i.e. the search estimator hugs its own
  min_ani=0.8 threshold on genuine ~90-95% ANI relatives. Smoke test: searching
  a genome against a DB containing **itself** returns ANI 0.9017 (4-enzyme DB) /
  0.9063 (BcgI DB) instead of ~1.0 — the GBRT v7 debias is grossly miscalibrated
  for near-identical inputs.
- Most of the 150 syn2bani hits are presumably spurious (cf. triangle result:
  0/2107 confirmed).

## 4. Accuracy cross-check: `dist` vs `ani` (500 stratified pairs)

Pairs drawn from `results/matrix_gtdb_r207_100k_v8_final.tsv`, stratified 125 per
band 80-85/85-90/90-95/95-100 by the validated v8 `ani` value
(`lists/accuracy_pairs.tsv`). Per pair: `ani` rerun, `dist --enzymes
BcgI,AlfI,AloI,FalI` (default chained-kmer output and `--mash-ani` GBRT v7),
default `min_af=0.1`. Full table: `dist_vs_ani.tsv`.

- **Sanity:** the `ani` rerun reproduces the matrix truth exactly
  (mean |d| = 0.0000 over 500 pairs) — the validated path is stable at 98177dc.
- **dist screens out 472/500 pairs (94.4%)** at default `min_af=0.1`: the coarse
  screen (exact shared tags / max tags) is miscalibrated for real draft-genome
  pairs whose true ANI is 80-100 but whose shared-exact-tag fraction is ~1-6%.
- On the 28 pairs that pass the screen (27 with ani>0): dist chained-kmer output
  is actually decent — Pearson r = 0.872, mean |dANI| = **0.85 points** vs the
  MLE truth (0.56 in the 95-100 band). The `--mash-ani` (GBRT v7) output is
  worse: mean |dANI| = 1.34 points, and it undershoots badly near identity
  (cf. the 0.90 self-hit in search).
- With `--min-af 0` the screen no longer filters, but the v7 matcher still
  returns ANI 0.0000 for some true 95-99% ANI draft pairs (3 of 4 spot-checked
  pairs, see `logs/smoke.out` debug section) — fragmented assemblies break the
  v7 chaining that `dist` relies on.

## Files

HPC `/lustre1/g/aos_shihuang/Syn2bANI-paper/results/db_scale/` and locally
`results/db_scale/` (tsvs + md + lists + logs only; sketch binaries stay on HPC):

- `sketch_scaling.tsv`, `triangle_scaling.tsv`, `search_scaling.tsv` — timing/RSS/store
- `dist_vs_ani.tsv` — 500-pair accuracy table (truth vs dist variants)
- `out/` — raw tool outputs (triangle edge lists, skani matrices, search results,
  `search_intersection_ani.tsv`, smoke tests)
- `lists/` — genome subset lists + accuracy pair list (reproducible, seed 42)
- `logs/` — `/usr/bin/time -v` logs, `cli_help.txt` (exact CLI surface),
  `slurm_3907622.out`, `smoke.out`
- `scripts/` — benchmark drivers + `analyze_search.py`, `analyze_accuracy.py`
- Sketch stores kept on HPC only: `sketches/s2ba_n5000` (491 MB),
  `sketches/skani_n5000` (1.7 GB)

## Bottom line for the manuscript

Only `sketch` can be benchmarked head-to-head with skani today, and it is
competitive (smaller store, linear scaling, ~8x slower wall-clock). `dist`,
`search`, and `triangle` need to be re-pointed at the validated
chain-restricted MLE estimator (or the v7 matcher needs to be fixed for
fragmented assemblies and the GBRT v7 head recalibrated near identity) before
any database-scale accuracy claims can be made. Specific repair list:

1. `triangle`/`search`/`db add`: stop hardcoding BcgI; honor the sketch's
   recorded panel (it is stored in `.s2ba`).
2. `dist`: add `--ql/--rl`; digest each genome once, not per query; recalibrate
   or replace the `min_af=0.1` exact-tag screen.
3. All three: replace the v7 `TagMatcher`/`AniCalculator` GBRT path with
   `chain_ani::compute`, or gate output behind a loud experimental warning.
4. `triangle` should emit NaN/empty, not `0.0000`, for below-detection pairs.
