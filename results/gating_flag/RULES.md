# Per-pair estimator gating and flag recalibration — rules and evidence

Date: 2026-08-14. Code: Syn2bANI `src/core/chain_ani.rs` (`gated_estimate`,
`unreliable`, constants `GATE_PARTIAL_GAP`, `FLAG_MAX_BP_PER_ANCHOR`) and
`src/cli/ani.rs` (new `ani_gated` / `gate` columns, new flag semantics).
All analysis scripts are in this directory; data joins in `prep_data.py`.

## Datasets

| set | n | truth | file(s) |
|---|---|---|---|
| GTDB-ANIm | 2,053 finite (of 2,074) | ANIm (dnadiff 1-to-1) | `gtdb_anim_joined.tsv` = `../anim_truth_2074_v8current.tsv` + `../anim_2074_acc2seqid.tsv` + `../panel_by_band/eval_pairs.tsv` |
| oral/gut | 100 same-species (`label==high`) | FastANI (reads ~1 pt low vs ANIm) | `oral_gut_joined.tsv` = `../oral_gut_1225_v8current.tsv` + mapping + `../../data/oral_gut_validation_merged_v8.tsv` |
| mid-ANI | 15 (truth 87.6–90.2) | ANIm | `midani_joined.tsv` + verbose 4-enzyme features rerun (`midani_15_verbose.tsv`) |
| sim ladder | 12 (85→99.9%) | exact by construction | `Syn2bANI/prototype/simindel*` |
| mosaic sims | 9 (90–98%) | exact by construction | `simulate_mosaic.py`, rerun locally (`/tmp/mosaic_test`) |

Likelihood-exact prototyping used per-enzyme strata dumped by the current
binary on HPC (`--strata-out`, 2,074 pairs, `strata_2074/`), with a Python
port of the Rust fits (`het_fit.py`) verified to reproduce `ani`,
`ani_uniform`, `retention` and the support/shape decision on a 40-pair
sample (`refit_cache.tsv` has LRT/NLL/shape for all pairs).

## Task A — gating rule

**Rule: report the homogeneous (uniform) estimate when
`|ani_from_loss − ani_from_hist| > 5` ANI points; otherwise the gamma
estimate (the existing `ani` column, including its LRT support gate).**
Shipped as the new `ani_gated` column; `gate` records
`gamma` / `uniform` / `uniform_fallback` / `none`.

Why this rule and not the alternatives (all evaluated in `task_a_oracle.py`
and the refit-cache analysis):

- **Oracle bound is tiny on GTDB.** Per-pair best-of-{gamma, uniform} gives
  MAE 2.788 vs always-gamma 2.881 — gating can win at most ~0.09 points
  there. The real payoff is the mid-ANI regime, where gamma overshoots low
  by 4–10 points.
- **BIC/LRT variants do not help.** BIC (`lrt > ln n_tags`): 2.913;
  chi-bar-square boundary threshold (`lrt > 2.706`): 2.881 (identical to
  current); stricter LRT (6/10/20): 2.92–3.05. The existing LRT support gate
  is already near-optimal on GTDB; the failure it misses is not a
  significance failure.
- **Significance-scaled discrepancy rules fail.** `|gamma−uniform| > k·SE`:
  3.46–3.50 — barely better than always-uniform (3.503). Retention gates and
  n_anchors gates: 2.89–3.43, all worse than always-gamma.
- **Effect-size gap wins.** `gap_lh ≤ τ` sweep: τ=3 → 3.161, τ=4 → 2.947,
  τ=4.5 → 2.862, τ=5 → 2.819, τ=5.5 → 2.811, τ=6 → 2.810, τ=7 → 2.826.
  Flat optimum 4.5–6; τ=5 chosen because mid-ANI needs τ ≤ 5.44 to catch
  12/15 pairs (mid-ANI gap range 3.60–6.36, median 5.70). Adding
  `retention < 0.2 ⇒ uniform` was tested and rejected (2.836 > 2.819:
  below-detection pairs with a small gap genuinely favour gamma, 1.39 vs 2.10).
- **The irreducible ambiguity:** on GTDB, pairs with retention 0.3–0.4 and
  gap 5–6 are *better* under gamma (by 1.26 MAE, n=16), while mid-ANI pairs
  in the same (retention, gap) cell are better under uniform by ~4.5 (n=12).
  No tested internal statistic (retention, het_shape/alpha, gamma-vs-loss
  position, n_anchors, af) separates these two groups — het_shape ranges
  overlap completely (mid-ANI 1.08–3.17 vs GTDB cell 1.52–4.28). τ=5 splits
  the difference; this is the honest limit of per-pair gating here.

Validation of the chosen rule (prototype, then confirmed with the new binary):

| set | always-gamma | always-uniform | gated | gate fires |
|---|---|---|---|---|
| GTDB-ANIm (MAE) | 2.881 | 3.503 | **2.819** (oracle 2.788) | 5.5% |
| mid-ANI (MAE) | 4.482 | 1.423 | **0.959** | 12/15 |
| oral/gut high (MAE) | 0.552 | 1.165 | **0.552** | 0/100 |
| sim ladder (MAE, exact truth) | 0.0732 | 0.0732 | **0.0732** | 0/12 |
| mosaic sims (MAE, exact truth) | 1.711 | 2.975 | **1.711** | 0/9 |

By band on GTDB (gated): 80–85: 3.369, 85–90: 3.490, 90–95: 1.817,
95–99: 1.159.

## Task B — flag rule

**New flag semantics** (`flag` column values unchanged):
`BELOW_DETECTION` = expected retention < 0.20 (unchanged, takes precedence);
`INCONSISTENT` = **gate fallback fired** (`|loss − hist| gap > 5` points)
**or `breakpoint_count / n_anchors > 0.5`**; else `ok`.

How it was chosen (`task_b_flag.py`, `task_b_gated.py`; AUC = ranking of
`|err_gated| > 1` on GTDB):

- Best single statistics: `bp_per_anchor` (AUC 0.803), `anchors_per_tag`
  (0.771), `std_err` (0.732). The old flag's mechanism — the loss-vs-hist
  gap at ~5 SE — has AUC 0.17–0.28 *inverted* on GTDB: any significance
  threshold of it flags the pairs the gamma correction helps. Regime
  conditioning (significance gap only at high retention) still inverts
  within every GTDB retention slice (e.g. retention ≥ 0.6: flagged 1.69 vs
  kept 3.61). That mechanism is unrecoverable as a reliability flag; it
  survives only as the effect-size gate/flag component above.
- The composite was checked against every dataset at the shipped thresholds;
  no inversion anywhere:

| set | kept MAE (n) | flagged MAE (n) |
|---|---|---|
| GTDB-ANIm | 2.415 (1,575) | 4.153 (478) |
| oral/gut high | 0.514 (99) | 4.269 (1) |
| mid-ANI | 0.343 (3) | 1.113 (12) |

  Old flag for contrast (gated errors): GTDB kept 3.883 / flagged 1.963
  (**inverted**), oral/gut 0.293 / 1.054, mid-ANI − / 0.959.

- **Honest cost:** on near-clonal oral/gut pairs the new flag is much less
  sensitive than the old one (kept-set MAE 0.514 vs 0.293) — the old flag's
  sensitivity there came entirely from the significance-gap mechanism that
  inverts on GTDB, and no threshold of it transfers. Likewise the mosaic sims
  (gamma bias up to +5.5) are now `ok` because neither component fires; the
  bias they exhibit is systematic, not per-pair detectable by these
  statistics.
- `retention < 0.35` as a flag component was tested and rejected: it inverts
  on GTDB post-gating (flagged 2.28 < kept 3.02 — high-retention GTDB pairs
  are the mosaic-biased ones).

## Task D — end-to-end before/after (new binary, HPC rerun)

New binary rerun over all 2,074 pairs on HPC (`../anim_truth_2074_gated.tsv`;
identical estimator columns to the previous matrix, plus `ani_gated`/`gate`
and the recalibrated `flag`). MAE vs ANIm (`before_after.tsv`):

| estimator | all | 80–85 | 85–90 | 90–95 | 95–99 | bias |
|---|---|---|---|---|---|---|
| gamma (before) | 2.881 | 3.571 | 3.529 | 1.817 | 1.159 | +2.41 |
| uniform | 3.503 | 4.098 | 4.271 | 2.381 | 1.568 | +3.50 |
| **gated (after)** | **2.819** | **3.369** | **3.490** | 1.817 | 1.159 | +2.68 |
| oracle bound | 2.788 | 3.274 | 3.478 | 1.802 | 1.150 | — |

Gate fallback fired on 112 pairs (5.5%), matching the prototype exactly.
New flag on the rerun matrix (gated error): ok n=1,531 MAE 2.445;
INCONSISTENT n=438 MAE 4.357; BELOW_DETECTION n=84 MAE 1.630. The ordering
is correct (flagged worse than kept). Note the BELOW_DETECTION set scores
*below* the ok set — those divergent pairs are where the gated estimator is
most helped by the fallback, and ANIm truth there is itself a partial-genome
measurement; the category means "no reliable estimate", not "large error".

Mid-ANI end-to-end with the new binary (`midani_15_gated.tsv`):
gamma 4.482 → gated **0.959** (uniform 1.423). The three pairs the gate keeps
on gamma are precisely the three where gamma was already correct
(errors 0.15 / 0.37 / 0.52); the other twelve fall back and are flagged
INCONSISTENT.

Sim ladder with the new binary: MAE 0.0732 (unchanged; gate never fires,
all rows `gate=uniform` because heterogeneity is correctly unsupported on
uniform-rate sims). Mosaic sims: MAE 1.711 (unchanged; gate keeps gamma on
all nine).

## Reproduction

```
python3 prep_data.py          # joins the three validation tables
python3 task_a_refit.py       # exact-likelihood refits from strata (needs strata_2074/)
python3 task_a_oracle.py      # oracle bound + simple gate sweep
python3 task_b_flag.py        # per-statistic AUCs (gamma error)
python3 task_b_gated.py       # per-statistic AUCs (gated error) + transfer
python3 task_d_eval.py        # before/after on the rerun matrix
```
