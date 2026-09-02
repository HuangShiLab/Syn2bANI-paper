# Syn2b structural channel on the SynTracker isolate cohorts

`08_syn2b_structural.py` runs Syn2b's structural channel (`breakpoints`,
`scj_distance`, `inverted_fraction`, `observable_fraction`) on the four isolate
collections from
[Enav, Paz & Ley, *Nat Biotechnol* 43:773–783 (2024)](https://doi.org/10.1038/s41587-024-02276-2)
and scores it against the answers that paper already published.

This is a different measurement from `04_run_syn2bani.sh`, which ran **Syn2bANI**
(ANI + `anchor_adjacency` + `breakpoint_count`). The structural channel has never
been run on these cohorts.

---

## Why these cohorts

Everything measured so far shows our estimator tracks dnadiff. Nothing shows it
stays **quiet** when there is nothing to find. These cohorts give us that, because
the paper pairs a synteny tool with a SNP tool and reports which mode dominates in
each collection:

| Cohort | n | Published result | What our channel must do |
|---|---|---|---|
| *Streptomyces rimosus* M527 | 20 isolates, 190 pairs | popANI pinned at 0.99990–1.0 (clonal) while APSS spans ~0.90–1.0 — variation is **structural** | **positive control**: signal must exceed the self-comparison floor |
| hypermutator *E. coli* | 23 isolates, 253 pairs | inStrain calls *no* pair the same strain; SynTracker calls *all* pairs the same strain — variation is **SNP-driven** | **negative control**: must stay at the floor |
| *Neisseria gonorrhoeae* | 12 isolates, 66 pairs | Spearman rho = 0.985 between the two scores — **both modes** | both channels move together |
| *Helicobacter pylori* | 77 isolates, 2,926 pairs | mixed; participants 322, 326, 439 carry subpopulations only one tool calls same-strain | the interesting case; connects to the cagPAI cohort |

The hypermutator *E. coli* set is the hard one: those genomes carry a heavy SNP
load, which is exactly what strips landmarks.

---

## Why the self-comparison control is mandatory

The Syn2bANI pass already in `data/syntracker_validation/syn2bani/` compared every
genome against **itself** as well as against the others. A genome is collinear with
itself, so those rows must read 0 breakpoints and adjacency 1.0. They do not:

| Cohort | self-comparison median `breakpoint_count` | cohort mean reported in `syntracker_summary.tsv` |
|---|---|---|
| *E. coli* hypermutator | 904.0 | 905.7 |
| *H. pylori* | 25.0 | 28.1 |
| *N. gonorrhoeae* | 327.5 | 305.3 |
| *S. rimosus* | 1415.0 | 1980.2 |

The reported per-cohort means reproduce each cohort's own self-comparison floor to
within a few percent. So the cross-species ordering in that table — read naturally
as "*S. rimosus* has the most structural variation" — is a ranking of **how
fragmented each cohort's assemblies are**. Three of the four cohorts carry
essentially no signal above the floor; *S. rimosus* shows some excess
(1980 vs 1415) but is still 71% floor.

The cause is the reference-side inflation in `breakpoint_count`, fixed in Syn2bANI
`c974f5f`. The point of this script is that **one control would have caught it
before the numbers reached a summary table**, so that control now runs first and
blocks the rest:

```
=== self-comparison control ===
cohort                              n   self bp  self SCJ  obs_frac    K_est  control
Test_cohort                         4       0.0       0.0    1.0000        1  PASS
```

A `FAIL` stops the run before any cohort statistic is written. `--allow-failed-controls`
overrides it for inspection, and stamps the failure into `self_control_<label>.tsv`.

**The old `results/syntracker_validation/syntracker_summary.tsv` and the old
`figures/syntracker_validation/` figures were derived from the pre-fix numbers and
have been replaced.** Updated summary, correlations and figures are now in place;
see `results/syntracker_validation/SYNTRACKER_STRUCTURAL_REANALYSIS.md`.

---

## What the metrics do and do not survive

Measured on a synthetic control built from *E. coli* K-12 — closed original; the
same sequence shattered into 40 contigs with half reverse-complemented; a 500 kb
inversion; a 1.2 Mb origin rotation — with the four-enzyme panel:

| Comparison | `breakpoints` | `scj_distance` | `inverted_fraction` | `observable_fraction` |
|---|---|---|---|---|
| closed vs rotated | **0** | **0** | 0.0000 | 1.0000 |
| closed vs 500 kb inversion | **2** | **4** | 0.1092 | 1.0000 |
| closed vs 40-contig shatter | **0** | 39 | 0.3535 | 0.9930 |
| shatter vs inversion | **2** | 43 | 0.3217 | 1.0000 |

Three things to take from this:

- **`breakpoints` is fragmentation-immune and rotation-invariant.** It reads 2 for
  every pair containing the inversion (2R for R = 1) and 0 for every pair without
  one, including the 40-contig shatter. This is the positive-contradiction rule: a
  contig break is an absence of evidence, not a contradiction.
- **`scj_distance` is not.** It carries a `+(K−1)` term by definition — the
  symmetric difference of adjacency sets genuinely loses an adjacency at each break
  (40-contig shatter → SCJ 39). **Do not read SCJ on draft assemblies.**
  `excess_over_floor` in the summary is therefore computed on breakpoints.
- **The orientation channel is only interpretable on near-closed assemblies.**
  `inverted_fraction` reads 0.3535 on the shatter, which is the drift toward 0.5
  predicted by the GTDB measurement (median `min(f, 1−f)` goes 0.2880 at K = 1–2 to
  0.4726 at K > 100). `observable_fraction` recovers K: 1 + (1 − 0.9930) × 5604 =
  40.2, against a true 40.

Truth for the inversion is 500,000 / 4,543,028 = 0.1101; measured 0.1092, a 0.8%
error.

---

## The reference-ordering step

These are short-read assemblies, so per the table above the orientation channel
needs help. The SynTracker paper hits the same wall and orders contigs against a
reference with Mauve before comparing.

Run the script **twice** and compare:

```bash
--label raw          # assemblies as produced
--label ref_ordered  # after Mauve contig ordering against samples/references.tsv
```

Reference-guided ordering biases each contig toward collinearity with the
reference, so inversions whose breakpoints fall on contig boundaries are absorbed
rather than detected. The bias is toward the null — it costs sensitivity, not
specificity — which is acceptable for a positive/negative control design but must be
declared in the paper rather than discovered in review. The gap between the two runs
is itself a measurement of what the orientation artifact costs.

---

## Running it

```bash
BASE=/lustre1/g/aos_shihuang/data/syntracker_validation
python3 $BASE/scripts/08_syn2b_structural.py \
    --assembly-dir $BASE/assemblies \
    --samples-dir  $BASE/samples \
    --out-dir      $BASE/syn2b_structural \
    --syn2b        /lustre1/g/aos_shihuang/Syn2b/target/release/syn2b \
    --label        raw \
    --workers      16
```

Roughly 3,435 within-cohort pairs plus 132 self-comparisons. Digestion is the
expensive step and is cached in `<out-dir>/tgt/`; a re-run with the same panel skips
it. Use a **separate `--out-dir` per enzyme panel** — the cache keys on isolate
alone.

### Outputs

| File | Contents |
|---|---|
| `self_control_<label>.tsv` | per-cohort self-comparison result, pass/fail, max deviation per metric, `K_est` |
| `syn2b_structural_pairs_<label>.tsv` | one row per pair, full structural channel, `same_group` flag for within-host pairs |
| `syn2b_structural_summary_<label>.tsv` | per-cohort medians, `self_floor_bp`, `excess_over_floor`, `K_est`, and the published expectation |

`excess_over_floor` is the only column that can carry biology. If it is ~0 for
*S. rimosus* the positive control has failed; if it is large for the hypermutator
*E. coli* set the negative control has failed.

---

## Two implementation notes

- **Text TGT, not binary.** `syn2b synteny` reads only the text path and does not
  sniff the format; a binary TGT fails with `stream did not contain valid UTF-8`.
- **The script renames each TGT's genome id to the isolate.** `syn2b digest` takes
  the genome id from the *first contig's FASTA header*
  (`main.rs`: `if gid.is_empty(){gid=rec.id.clone();}`). On SPAdes assemblies that
  is `NODE_1_length_..._cov_...`, which makes `genome_A`/`genome_B` unreadable and,
  if two isolates ever share a first-contig name, makes `syn2b synteny` refuse the
  pair outright. The same rename is what makes the self-comparison possible at all,
  since Syn2b rejects two genomes sharing an id.

## Prerequisites

`samples/samples_<cohort>.tsv` (columns `species`, `isolate`, optionally
`host`/`participant`, `sra_run`) and `assemblies/<isolate>.fna`, both produced by
steps 00–03 of the main [README](README.md).
