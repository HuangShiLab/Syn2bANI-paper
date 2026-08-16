# IUPAC geometry fix — validation report (2026-08-17)

Status: IMPLEMENTED and deployed (code repo commit f054dbb). Companion to
DESIGN.md, which it follows; one forced deviation, documented below.

## Deviation from the design premise: digestion was not correct

DESIGN.md assumed "digestion is already correct; only the MLE geometry is
wrong." A direct probe disproved this: the static digest patterns for the
three degenerate enzymes enforced only the fixed *prefix* of each anchor plus
the degenerate bases, dropping the trailing fixed bases of the
degenerate-containing right anchor:

- HaeIV `GAY-N5-RTC`: only `GA`, `Y`, `R` enforced; trailing `T`,`C`
  unconstrained (breaking TC still emitted a tag; breaking Y removed it)
- Hin4I `GAY-N5-VTC` (rev `GAB-N5-RTC`): trailing `TC` unconstrained on both
  strands
- BaeI `AC-N4-GTAYC` (rev `GRTAC-N4-GT`): final `C` (fwd) / `GT` (rev)
  unconstrained

Consequence: those positions behaved as body in reality, so even the perfect
DESIGN.md geometry would have stayed biased against the actual digest
(measured +0.29% at 95% ANI on a digest-level simulation with geometry fixed
but digest permissive; old geometry +0.35% — the geometry fix alone would
have looked like a no-op). Historical data corroborates the permissive
digest: HaeIV yielded 13,319 mean shared tags vs BcgI's 316 (density
~1/64/strand, not the full site's 1/1024).

Fix (src/enzyme/digest.rs): added the missing fixed anchors — HaeIV `TC`,
Hin4I `TC` (both strands), BaeI `C` (fwd) / `GT` (rev). Digestion now enforces
the complete biological sites, matching the registry anchor strings and the
MLE geometry. The default panel (BcgI/AlfI/AloI/FalI) patterns are untouched.

## As implemented

- `SiteGeometry { tag_len, exact_site, d2, d3 }` per enzyme, derived by
  parsing the IUPAC anchor strings (A/C/G/T exact; Y R S W K M → d2;
  V H D B → d3; N → body). Sketch-side fallbacks route through
  `geometry_for_name` (unknown enzyme → historical 32/6/0/0 default).
- Homogeneous likelihood: the exact convolution of DESIGN.md, accumulated by
  log-sum-exp. Heterogeneous: NB with effective body `b + d2/3 + 2·d3/3` plus
  the closed-form survival factor `α·ln(α/(α + d·(2·d2/3 + d3)))`.
- At d2 = d3 = 0 both paths are bit-identical to the pre-fix code (guarded
  branch + a `to_bits` unit test), so every result with the default panel is
  unchanged.
- BaeI's exact-site count is 6 (AC + GTAC); DESIGN.md's "exact 5" was an
  arithmetic slip in the doc text.

## Validation

- `cargo test --release`: 96 green (85 lib + 11 integration). New tests:
  - Monte-Carlo count recovery at d2/d3 > 0, simulation independent of the
    likelihood code: homogeneous within 2e-4 of truth at ANI 0.85–0.98
    (HaeIV-like e4/d2=2, Hin4I-like e4/d2=1/d3=1); the same counts under the
    old geometry give +0.19pp bias at 0.90, +0.22pp at 0.85. Het-model counts
    at (d, α) ∈ {(0.02, 0.5), (0.05, 1.0), (0.10, 2.0)} recovered within 2e-3.
  - Digest-level: 2 Mb genome, 5% substitutions, real HaeIV digestion +
    chaining: new geometry 0.95175 vs truth 0.95 (residual ≈ 1 s.e.
    Monte-Carlo noise).
  - Digest regression: every exact-site mutation kills the tag, within-class
    degenerate mutations preserve it (all three enzymes, both strands).
- simindel regression (12 pairs, default panel): **byte-identical
  before/after** — the regression gate holds exactly.
- simindel with `-e HaeIV,Hin4I` (truth vs estimate, ANI%):

  | truth | before (old digest+geometry) | after |
  |---|---|---|
  | 85.0 | 85.62 (+0.62) | 84.37 (−0.63) |
  | 90.0 | 90.63 (+0.63) | 89.67 (−0.33) |
  | 95.0 | 95.28 (+0.28) | 94.80 (−0.20) |
  | 98.0 | 98.03 (+0.03) | 97.99 (−0.01) |
  | 99.9 | 99.89 (−0.01) | 99.90 (+0.01) |

  MAE 0.273 → 0.242 pp. The upward geometry bias is gone; the small residual
  *downward* drift at high divergence is consistent with the simindel indel
  fraction (kept_frac ≈ 0.989 → tag loss the mutation-only model reads as
  divergence), not with site geometry.

## Caveats / follow-ups

- `.s2ba` sketches and `panel` strata files containing HaeIV/Hin4I/BaeI
  produced before this change reflect the old permissive digest (~16× more
  tags) and must be rebuilt before reuse with the corrected geometry.
- Any paper-repo result measuring HaeIV/Hin4I/BaeI tag density or panel
  sweeps including them (e.g. the 11-enzyme panel MAE 0.670, benchmark
  timings) predates both fixes; re-run the enzyme-panel sweep at high
  divergence before citing.
- The digest fix cuts HaeIV/Hin4I tag density ~16×; their standalone
  estimator power drops accordingly. Density is a panel-design question, not
  a correctness one.
- Unblocks: the wide-panel revisit and the capped-NPMLE revisit
  (fits cached in results/spatial_model/gtdb_spatial.tsv).
