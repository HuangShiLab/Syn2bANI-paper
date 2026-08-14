# Spatial/mechanistic models for the divergent-pair ANI bias — candidates, evidence, verdict

Date: 2026-08-14. Code: Syn2bANI `22b853b` (unchanged — this is a Python-only
Phase 1). Prototype scripts and per-pair result TSVs live in this directory.
Primary data: `../gating_flag/gtdb_anim_joined.tsv` (2,053 finite GTDB-ANIm
pairs), `../gating_flag/strata_2074/` (per-enzyme sufficient statistics),
19 exact-truth mosaic sims (9 canonical from `simulate_mosaic.py` + 10 extra
replicates with other seeds, block sizes 1/5/20 kb, ANI 0.85–0.95),
12 uniform-rate sims (`simindel`), mid-ANI 15 pairs (ANIm), oral/gut 100
same-species pairs (FastANI).

**Verdict up front: negative result.** No mechanistic candidate beats the
gated baseline (GTDB MAE 2.819) meaningfully while holding the simulation
gates, and none comes close to the ridge calibration (0.988). The bias
information the ridge uses is **statistical, not mechanistic, at this tag
density**: the dominant error term lives in the unchained genome fraction,
whose divergence is not identifiable from tag data. The mechanism of the bias
is nonetheless now established quantitatively, and one sub-component (the
in-chain distributional distortion) is genuinely fixable — it is just worth
only ~0.1–0.2 MAE on GTDB.

## 1. Mechanism, established on exact-truth sims

Running the current binary on the 19 mosaic sims and comparing chain coverage
(PAF spans) against the known per-block rates (`explore_mosaic.py`,
`measure_pi.py`):

- Chain coverage of a 5 kb block declines softly with its divergence v:
  ~0.99 at v<0.07, ~0.8–0.95 at v≈0.1–0.2, ~0.6 at v≈0.3. It is **not** a
  step: moderately divergent blocks are carried into chains by their anchored
  flanks (the pass-2 skip budget at 90% ANI tolerates ~9 consecutive tag
  failures ≈ 40 kb). Coverage is strongly scale-dependent: at 1 kb blocks AF
  ≈ 0.99 even for bimodal 90% ANI, at 20 kb blocks AF ≈ 0.76 for the same
  rate distribution.
- The bias decomposes into two additive terms:
  1. **Coverage term.** The chained fraction is more conserved than the whole
     genome. An oracle that knew the exact mean identity of the chained
     sample still errs by MAE 1.359 (bias +1.359) against whole-genome truth
     on the 19 sims.
  2. **Family term.** Even within the chained sample, the parametric fits
     misplace mass: on *asymptotic* in-chain counts (expected counts from the
     true per-tag divergences, no sampling noise; `asymptotic_floor.py`),
     gamma errs by MAE 2.251 (+2.164), a free two-component mixture by 2.204,
     a capped grid NPMLE by 1.252. The 4 categories per enzyme (m = 0, 1, 2,
     miss) leave the tail mass weakly identified: the fit trades mid-range
     mass against deep-tail mass almost freely, and the gamma family resolves
     that trade wrongly under ascertainment.
- The b1k replicates isolate the family term: coverage ~complete (AF 0.992)
  yet gamma reads +1.1 (gamma-rates case) to +4.2 (bimodal case) high.

## 2. The unchained fraction is not identifiable from tag data

The only genome-wide observable that does not share the chain denominator is
the anchor count. On the worst mosaic case (bimodal 70% core, 90% ANI):

- `n_anchors` = 4,144; found-in-chains = 3,567; the excess (577) is fully
  explained by multi-match inflation κ ≈ 1.16 (repeat/duplicate matches),
  leaving an out-of-chain match residual of ≈ 0 — while the true unchained
  tags (829) have match probability ≈ 0.02 and contribute ≈ 7–20 expected
  anchors. **The signal (~0.2% of tags) is below the multi-match noise**, so
  the anchor residual can only say "the unchained mass is saturated-divergent
  *or* accessory" — the two cases are indistinguishable.
- Consequence for candidate (a): `ani_unchained` cannot be estimated per pair
  from the loss side. AF-weighted mixtures must supply it from an assumption.

On real GTDB pairs the unchained fraction is accessory-dominated, and ANIm
excludes accessory too. The required "unchained identity"
`a_u = (anim − AF·gated)/(1−AF)` is tightly distributed but **band-determined,
not pair-determined**: medians 84.9 / 86.9 / 90.9 in the 80–85 / 85–90 /
90–95 bands (p10–p90 ≈ ±1.7–3), i.e. only ~0.5–1.5 points below the band mean
truth — nothing like the 65–77% divergent-tail identity of the mosaic sims.
Within every band, truth is uncorrelated with AF (r ≈ +0.06…+0.20) while the
estimator's error correlates with AF at r ≈ −0.72…−0.80 (slope +4.5 to +8.4
ANI points per unit unchained fraction). That association is exactly what the
ridge calibration learns; it is not recoverable from per-pair tag statistics.

## 3. Candidates and results

Derivations and fits: `model_spatial.py` (C1–C3, A1), `model_discrete.py`
(D2–D4, NPcap), `gate_lrt.py` (LRT gate), `eval_gtdb.py`, `eval_controls.py`.

- **(a) AF-weighted mixtures.** Several algebraically distinct variants were
  written down and tested: identity space `AF·ani_c + (1−AF)·a_u`, divergence
  space `−ln` of the same, and match-probability space (combine E[q], invert
  once). With `a_u` from the anchor residual (`A1`): mosaic MAE 1.527 —
  best closed form — but unstable across heterogeneity scale (−6.4 at 20 kb
  blocks where `a_u` hits its detection floor, +4.2 at 1 kb). On GTDB the
  mosaic-calibrated `a_u = 0.70` is catastrophic (MAE 11.1): the unchained
  mass there is accessory. `a_u = 0.85` *looks* good on GTDB (MAE 1.964) but
  fails the mosaic gate (+4.1 on the 1 kb bimodal case) — it is fitting the
  ANIm partial-coverage artifact, not mechanism. **Rejected by the gate.**
- **(b) Discrete two-component / K-component ML.** D2: mosaic 1.708 (+1.07),
  but uniform-sim MAE 0.895 — invents heterogeneity on uniform data,
  **fails the ≤0.1 gate outright**. D3/D4 are erratic (single-case blowups to
  −10). BIC prefers D2 over gamma on only 5/19 mosaic cases — per-pair data
  do not support the extra component. The ascertainment-aware parametric
  fits (C2: 2-component + coverage tilt + AF constraint, 2.757; C3: gamma +
  coverage tilt + AF constraint, 2.554) **overshoot to −5 on several cases** —
  a third instance of the tilted-Gamma failure pattern. The coverage function
  π(v) is soft and context-dependent (§1), and these fits are acutely
  sensitive to it. Not tunable honestly; abandoned.
- **(d) Capped grid NPMLE + LRT gate** (the strongest candidate). A 66-point
  divergence grid over v ∈ [0, 0.45] — the cap is the principled statement
  that beyond ~35% divergence a homologous region is indistinguishable from
  accessory at 32 bp tags — fitted by EM to the in-chain counts; ANI = mean
  identity under the fitted distribution. Gate: use NPMLE only when
  `2·(NLL_gamma − NLL_npmle) > 5`, else the gated estimate.

### Per-dataset results (MAE / bias, ANI points)

| estimator | mosaic 19 (exact) | uniform 12 (exact) | GTDB 2,053 (ANIm) | mid-ANI 15 (ANIm) | oral/gut 100 (FastANI) |
|---|---|---|---|---|---|
| gated (baseline) | 2.642 / +2.57 | **0.073** / +0.06 | **2.819** / +2.68 | **0.959** / +0.89 | 0.552 / +0.42 |
| NPMLE capped | **1.474** / +0.16 | 0.216 / −0.19 ✗gate | 2.707 / +2.04 | 1.525 / −1.02 | **0.425** / +0.34 |
| NPMLE, LRT-gated | 1.793 / +0.57 | **0.058** / +0.05 | 2.867 / +2.61 | 1.359 / −0.86 | 0.535 / +0.38 |
| 50/50 average | 1.725 / +1.37 | 0.107 / −0.06 | 2.634 / +2.36 | ~interpolates | ~interpolates |
| D2 (two-component) | 1.708 / +1.07 | 0.895 / −0.88 ✗gate | — | — | — |
| ridge calibration | — | — | 0.988 | 1.34 (Set A) | 0.460 (Set A) |

GTDB by band (gated / NPMLE / LRT-gated / average): 80–85: 3.369 / 3.143 /
3.384 / 2.978; 85–90: 3.490 / 3.346 / 3.528 / 3.314; 90–95: 1.817 / 1.802 /
1.892 / 1.741; 95–99: 1.159 / 1.220 / 1.257 / 1.150.

Reading: the capped NPMLE genuinely fixes the family term (mosaic 2.64 →
1.47, close to the 1.25 asymptotic floor) and helps slightly everywhere on
real data, but it fails the uniform-sim gate ungated (0.216 > 0.1) and
degrades mid-ANI (1.525 vs 0.959 — the LRT fires on the same low-retention
structure the effect-size gate already handles better). With the LRT gate it
passes uniform sims (0.058) but then *loses* to the gated baseline on GTDB
(2.867 vs 2.819). The plain average is the best GTDB variant (2.634) — a
0.185 gain that does not survive the gates either (uniform 0.107, mid-ANI
interpolates toward 1.2). Nothing satisfies "beats 2.819 meaningfully AND
holds every gate".

## 4. Why per-block dumps (roadmap option c) are not warranted

Per-block match/miss counts (`ani --blocks-out`) would sharpen the *in-chain*
rate distribution — the family term, worth ≤0.6 MAE on GTDB (2.819 → ~2.7
best case, the gap between the gated fit and the chained-sample oracle). They
cannot address the coverage term, which lives outside the chains and is where
the GTDB bias concentrates (§2). The skip: no Rust addition. If tag density
ever rises (larger enzyme panels), the calculus changes: deeper retention
makes the divergent tail directly observable in-chain, and the capped-NPMLE
family fix is the first thing to revisit.

## 5. Standing conclusions

1. The +2.4–2.7 GTDB bias is ascertainment, mechanistically confirmed: within
   every ANI band the error tracks AF (r ≈ −0.75) while truth does not.
2. Its pair-level magnitude is **not** encoded in any tag statistic that
   exists at 4-enzyme density: the unchained mass's divergence is observed
   only through an anchor residual that is zero within multi-match noise, for
   saturated-divergent and accessory DNA alike. ANIm's own partial coverage
   (48–72% of genome at divergent bands) makes the truth itself a
   conserved-weighted measurement, so part of the "bias" is an estimand
   mismatch no internal statistic can close.
3. The ridge calibration's 0.988 therefore does not have a mechanistic
   competitor at this tag density. This is the documented negative result the
   roadmap allowed for: the bias information is statistical (band/AF-level),
   not mechanistic (pair-level).
4. One real, transferable improvement exists — the capped-NPMLE in-chain fit
   removes the family distortion (mosaic 2.64→1.47; GTDB −0.11; oral/gut
   −0.13) — but it fails the uniform-sim gate ungated and the mid-ANI control
   gated, so it is **not shipped**. It is recorded here with its fits
   (`gtdb_spatial.tsv`) for reuse if tag density increases.

## 6. Reproduction

```
# sims (regenerates genomes; needs Syn2bANI checkout + release binary)
python3 /Users/macstudio/Downloads/Syn2bANI/prototype/simulate_mosaic.py \
    /Users/macstudio/Downloads/Syn2bANI/prototype/mg1655.fasta mosaic
python3 measure_pi.py --gen-extra      # extra replicates + coverage function
python3 explore_mosaic.py              # mechanism, section 1
python3 asymptotic_floor.py            # infinite-data floors
python3 model_spatial.py               # candidates C1-C3, A1
python3 model_discrete.py              # D2-D4, capped NPMLE
python3 gate_uniform.py                # uniform-sim gate
python3 gate_lrt.py                    # LRT gate landscape
python3 eval_gtdb.py                   # GTDB-ANIm (writes gtdb_spatial.tsv)
python3 decompose_gtdb.py              # section 2 diagnostics
python3 eval_controls.py               # mid-ANI + oral/gut controls
python3 final_summary.py               # summary_table.tsv, bic_mosaic.tsv
```

Mid-ANI strata were produced on HPC (binary at `32294d3`, strata identical to
`22b853b` — the intervening commits touch gating columns and the breakpoint
formula only) from `/lustre1/g/aos_shihuang/data/validation_mid_ani/genomes/`;
oral/gut strata are the lab's `results/oral_gut_strata.tsv`. Mosaic/simindel
query FASTAs are regenerable and not kept in the repo; all estimator outputs
(`.ani.tsv`, `.paf`, strata) are kept.
