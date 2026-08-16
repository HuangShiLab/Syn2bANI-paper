# IUPAC geometry fix — design note (2026-08-15)

Status: DESIGNED, not implemented. Implement after the db-path rewrite lands
(agent-12 owns the code repo until then).

## Problem

`geometry_from` (src/core/chain_ani.rs:206) computes
`site = left_anchor.len() + right_anchor.len()`, treating every IUPAC
degenerate position as fully specific. Digestion itself is correct
(`digest.rs` enforces `IupacConstraint` bitmasks), so tag *existence* is
right; only the MLE geometry is wrong. Affected enzymes:

- HaeIV  GAY / RTC  -> exact 4 (GA,TC), degenerate 2-of-4 x2 (Y, R); current site_len 6
- Hin4I  GAY / VTC  -> exact 4 (GA,TC), 2-of-4 x1 (Y), 3-of-4 x1 (V); current site_len 6
- BaeI   AC / GTAYC -> exact 5 (AC,GTAC... check: G,T,A?,C — Y only degenerate), 2-of-4 x1 (Y); current site_len 7

A mutation inside a degenerate class (C<->T at a Y position) preserves the
recognition site: the tag survives and the position shows up as a *mismatch*
in the histogram (tags are compared as literal 2-bit bases). The current
model assigns these positions survival `a` and zero mismatch rate, so it
overstates the site constraint and understates the body. This is what blocks
the wide panel at high divergence and gates the capped-NPMLE revisit
(results/spatial_model/MODELS.md).

## Correct model

Per position type (identity framework, mutant base uniform over 3
alternatives):

- exact site (e positions): survival a, never a mismatch
- body (b positions): survival 1, mismatch w.p. (1-a)
- degenerate k-of-4 (d2 = count of Y/R/S/W, d3 = count of V/H/D):
  survival a + (k-1)(1-a)/3; site-preserving mutation w.p. (k-1)(1-a)/3
  (counted as a mismatch); site-killing mutation w.p. (4-k)(1-a)/3
  (tag gone)

Homogeneous likelihood (exact convolution, d2+d3 <= 3 so cheap):

  P_m(a) = a^e * sum over j2,j3 with j2+j3<=m of
    C(d2,j2) q2^j2 a^(d2-j2) * C(d3,j3) q3^j3 a^(d3-j3)
    * C(b, m-j2-j3) (1-a)^(m-j2-j3) a^(b-(m-j2-j3))
  where q2 = (1-a)/3, q3 = 2(1-a)/3.

Heterogeneous (gamma-mixed) likelihood: given regional rate r, independent
Poisson channels — body mismatches rate r*d*b, preserving mutations rate
r*d*(d2/3 + 2*d3/3), killing mutations rate r*d*(2*d2/3 + d3/3). Gamma
mixing keeps them independent, so:

  ln P(m) = ln NB_pmf(m; alpha, eff_len = b + d2/3 + 2*d3/3)
          + alpha * ln( alpha / (alpha + d * (2*d2/3 + d3) ) )

i.e. the degenerate positions extend the NB effective body length and add a
closed-form survival factor. Reduces to the current model when d2=d3=0.

## Implementation plan

1. `Geometry` value becomes a small struct { tag_len, exact_site, d2, d3 }
   (or keep the tuple and add a side table keyed by enzyme name — check
   call sites: chain_ani.rs:1003, cli/ani.rs geometry_from_sketches,
   cli/struct.rs:123).
2. Derive the degenerate counts from the enzyme's IUPAC anchors in
   `geometry_from` (parse left/right anchor strings: Y R S W -> d2,
   V H D -> d3, N -> body).
3. mle.rs: replace the binomial/`a^site` terms with the convolution and the
   NB + survival-factor forms above, in both the homogeneous and het paths
   (and the per-enzyme agreement fit).
4. Tests: synth-count recovery at d2/d3 > 0 (homogeneous + het), plus a
   digest-level simulation on a pair built with HaeIV/Hin4I tags.
5. Validation gate: uniform sims (must stay ~unbiased), mosaic sims,
   mid-ANI, GTDB-ANIm 2,074 pairs — no regression elsewhere; then re-run
   the enzyme-panel sweep at high divergence and revisit capped NPMLE
   (fits cached in results/spatial_model/gtdb_spatial.tsv).

Watch out: `mismatches()` caps comparisons at 32 bases (2-bit packing);
stratum() already mirrors this with tag_len.min(32) — keep that interplay.
