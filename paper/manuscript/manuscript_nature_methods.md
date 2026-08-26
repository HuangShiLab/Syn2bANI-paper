# Syn2b-ANI: strain-level ANI and structural comparison from fixed restriction-site anchors

## Authors

Yufeng Zhang^1^, Shi Huang^2^

**Affiliations**
1,2 Faculty of Dentistry, the University of Hong Kong, Hong Kong SAR, China

**Corresponding author**
Shi Huang — shihuang@hku.hk

---

## Abstract

Average Nucleotide Identity (ANI) is the standard metric for prokaryotic species delineation, yet alignment-free tools reduce a genome comparison to a single divergence scalar and discard synteny and structural information. We introduce Syn2bANI, which uses in-silico Type IIB restriction enzyme tags as fixed positional anchors, chains shared anchors collinearly, and estimates ANI by maximum likelihood from tags inside trusted chains. A per-pair effect-size gate chooses between uniform and gamma rate models, and a ridge calibration trained on alignment-based truth corrects real-genome rate heterogeneity. In one pass the same anchors yield aligned fraction, a synteny score, and structural-variant calls. On simulations with exact truth Syn2bANI attains MAE 0.073% versus 0.377% (skani) and 0.742% (FastANI), and remains robust to indels, GC variation, accessory content, and fragmentation. On 43,334 held-out GTDB-R207 pairs it achieves MAE 0.619% against ANIm, compared with 0.957% for skani and 0.977% for FastANI. Syn2bANI is implemented in Rust and available at https://github.com/HuangShiLab/Syn2bANI.

**Keywords**: ANI, restriction-site anchors, structural variation, synteny, maximum likelihood, metagenome-assembled genomes

---

## Introduction

Average Nucleotide Identity (ANI) — the mean nucleotide identity of orthologous regions shared by two genomes — underpins prokaryotic species delineation and strain-level relatedness1,2. It is computed at scale in dereplication pipelines, outbreak investigations, and metagenomic surveys, increasingly on metagenome-assembled genomes (MAGs) that are fragmented by construction3.

Three tool families dominate. Alignment-based methods (ANIm/ANIb) are the direct realization of the definition but too slow for database-scale all-vs-all screens4,5. FastANI popularized an alignment-free approximation from 20-mer seed matches and fragment-level chaining2; skani replaced seeding with sparse k-mer sketches and a regressed chain statistic, gaining orders of magnitude in speed6. Both k-mer tools reduce the comparison to essentially one number: ANI plus an aligned-fraction estimate. Neither reports *where* the genomes agree — collinear blocks, inversions, translocations, or indels — even though that structural information is computed and then discarded inside their own chaining steps.

Type IIB restriction enzymes offer a different anchoring substrate. They cleave at fixed offsets from a short recognition sequence, producing short (~32 bp) tags whose genomic positions are fully determined by the motif. The 2bRAD reduced-representation sequencing family exploits this property7, and 2bRAD-M has shown species- and strain-level resolution in real microbiome data8. Computationally, these tags are sparse (one per ~0.5–1.5 kb), long enough that a match is almost certainly homologous, and — crucially — the *absence* of a tag at a homologous position is informative, implying a mutation in the recognition site or tag body.

Here we present Syn2bANI, which builds ANI estimation on these fixed anchors. Both genomes are digested in silico with a panel of Type IIB enzymes (default BcgI, AlfI, AloI, FalI). Shared tags become anchors; anchors are chained by a gap-penalized collinear dynamic program; and ANI is estimated by maximum likelihood from per-enzyme mismatch and miss counts *inside chains only*. Restricting the likelihood to chained regions separates divergence of the shared genome (ANI) from the fraction shared (aligned fraction). Because anchors carry coordinates, the same pass yields synteny scores and structural-variant calls at no extra cost.

We validate Syn2bANI in two regimes. Simulations with exact truth establish that the raw estimator is essentially unbiased and quantify robustness to indels, fragmentation, accessory content, and GC. Real genome pairs against ANIm truth reveal regional rate heterogeneity that biases the raw estimator upward by about two ANI points at divergent bands; a transparent ridge calibration on internal features removes this bias. The result is a tool that matches or exceeds skani and FastANI on held-out GTDB-R207 pairs while additionally providing synteny and SV outputs that k-mer tools cannot.

---

## Results

### The Syn2bANI estimator

Syn2bANI compares two genomes in four stages (Fig. 1a). First, both genomes are digested in silico with the default four-enzyme panel (BcgI, AlfI, AloI, FalI). Second, tags shared between genomes — exact or within a mismatch budget of 2 — become anchors; pigeonhole seeding makes tolerant recovery cheap. Third, anchors are grouped by (query contig, reference contig, orientation) and chained by gap-penalized collinear dynamic programming in two passes, with the second pass using a chain-break threshold derived from a provisional divergence fit. Fourth, per-enzyme tag outcomes inside chains enter a maximum-likelihood estimator of per-base divergence.

The likelihood treats a Type IIB tag as a tag that contains its own recognition site: a mutation in the site deletes the tag (miss), while mutations in the tag body appear as mismatches up to the budget. The mismatch histogram and miss rate are combined by Fisher information. Two nested models are fitted — a uniform single-rate model and a gamma rate-heterogeneity model — and an effect-size gate chooses per pair. Two diagnostics accompany every estimate: a consistency flag (`ok` / `INCONSISTENT` / `BELOW_DETECTION`) and a one-sided 95% upper bound for pairs below the detection floor. Full details are in the Online Methods.

The design decouples ANI from aligned fraction by construction. Accessory content outside chains is excluded from the ANI denominator, while AF reports the shared fraction independently. A short-contig rescue pass ensures that contigs carrying ≥2 placeable unique anchors but no accepted chain still contribute to the likelihood (Online Methods), which is essential for real draft assemblies.

### Accuracy and robustness under exact truth

We first validated the estimator on simulations in which true ANI is known exactly: query genomes are evolved from *Escherichia coli* K-12 MG1655 by applying a counted number of substitutions, and every ladder genome carries a 400 kb inversion plus ~46 deletions of 200–2,000 bp. Over twelve pairs spanning 85.0–99.9% true ANI, Syn2bANI attains MAE **0.073** ANI points, against **0.377** for skani and **0.742** for FastANI (Fig. 1b). Both comparators underestimate systematically toward low ANI. The 85% rung is flagged `BELOW_DETECTION` rather than reported as confident, by design.

Four further exact-truth families isolate confounders (Fig. 2; Supplementary Figs. S5–S8). At fixed 95% ANI with 0–4 deletions per 100 kb, errors stay between +0.032 and +0.184 (MAE 0.081). Simulated fragmentation into 20–201 contigs yields errors of +0.018 to +0.179 (MAE 0.093). Ladders on five genomes spanning GC 27.2–72.1% give MAE 0.074–0.356, not monotone in GC. With 0–50% of the genome replaced by homology-destroying shuffled blocks, the estimate stays flat (MAE 0.114) while `af_query` tracks the true shared fraction to within 0.004. Real *E. coli* draft assemblies confirm the short-contig rescue: bias at N50 5 kb improves from +0.65 points before the fix to +0.20 after, with AF recovering from 0.434 to 0.580.

### Independent-truth benchmark against GTDB-R207

The main accuracy result uses GTDB-R207 genome pairs with ANIm (dnadiff 1-to-1) as truth. A 2,074-pair training set across ANI bands (80–85%, 85–90%, 90–95%, 95–99%) was used for calibration; a separate 467-pair expansion at 95.0–99.15% was added because the original top band was underpowered.

On real genomes the raw estimator overestimates divergent pairs: post-rescue raw gamma MAE is 2.35 against ANIm, worst in the 80–85% band (MAE 3.21). Mosaic simulations with exact truth show the mechanism: real pairs mix a conserved core with more divergent regions, tags survive preferentially in conserved regions, and any single-rate fit reads high. The gamma model captures part of this bias, but a residual remains. A ridge regression on nine internal features (gated ANI, uniform ANI, both AFs, standard error, retention, anchor/chain/tag counts) removes the remainder. The model is trained by band-holdout cross-validation: each ANI band is predicted by a model trained only on the other three bands.

The deployed calibration v5 attains band-holdout CV MAE **0.731** (bias +0.05, r = 0.963) on the 2,074-pair benchmark, ahead of skani (0.906) and FastANI (1.056) (Fig. 3a). On a strictly held-out set of **43,334** same-genus pairs with training genomes excluded in both directions, the calibrated estimator attains MAE **0.619** (bias −0.12, r = 0.962), versus skani 0.957 and FastANI 0.977 on the same pairs (Fig. 3b). The advantage is largest at 80–90% ANI, where skani and FastANI carry a ~1.9-point downward bias. The one-sided upper bound covers ANIm truth on 99.4% of all 43,334 pairs.

A unified 80–100% benchmark combines the 43,334 held-out pairs with 727 high-ANI test pairs sampled from non-representative GTDB-R207 genomes. A simple hybrid rule — raw gated estimate at ≥98% ANI, calibrated estimate below — gives the lowest overall MAE (**0.615**), recovering the strong raw performance in the 97–100% regime (MAE 0.23) without degrading 95–97% (Fig. 4; Supplementary Fig. S3).

### Accuracy on metagenome-assembled genomes

On 695 CAMI2 MAG bins paired with their dominant source genomes, the raw gated estimate attains MAE **0.163** against dnadiff ANIm, versus 0.092 for skani and 0.203 for FastANI (Fig. 5). By CheckM2 tier, MAE is 0.041 on high-quality bins (n = 200; 100% within 0.5 points), 0.105 on medium-quality, and 0.291 on low-quality; |error| tracks contig N50 at Spearman ρ = −0.68. Applying the GTDB-trained v5 calibration to these MAGs degrades accuracy sharply (MAE 1.26), because the calibration learned real-genome associations absent in draft inputs; the raw gated estimate is therefore recommended for MAGs.

### High-ANI pairs hide extensive rearrangements

ANI is a scalar average and therefore blind to genome architecture. We quantified this on four published isolate collections. In a longitudinal *E. coli* hypermutator lineage, every cross-time pair with ANI >99.9% shows a synteny score of only 0.70–0.76 and 850–1,130 breakpoints. Re-analysis with `syn2bani struct` resolves representative pairs into dozens of chains with tens of inversions and translocations (Fig. 6a; Table 1). *Helicobacter pylori* same-host pairs cluster at ANI >99.9% but the lowest-synteny cases carry 68–103 breakpoints. *Neisseria gonorrhoeae* and *Streptomyces rimosus* near-clonal pairs show the same decoupling, with *S. rimosus* reaching 2,957 breakpoints at 99.9945% ANI (Fig. 6b; Supplementary Fig. S11).

The same pattern appears at database scale. Among 336 GTDB-R207 representative pairs with ANI ≥ 99.0% and aligned fraction ≥ 0.9 in both directions, 7.4% have a synteny score < 0.99 and 2.4% fall below 0.98 — i.e. pairs an ANI-only search would return as near-clonal yet which are structurally divergent.

### Computational efficiency

On all-vs-all comparison of 22 genomes (14 complete Enterobacteriaceae plus 8 real draft *E. coli*), Syn2bANI computes all 484 pairs in 3.73 s from FASTA or 2.79 s reusing sketches, against 17.7 s for FastANI and 0.15 s for skani (the latter reports only 302/484 pairs because of its minimum-AF filter and excludes sketching). Peak memory is 58 MB with sketch reuse versus 185 MB (skani) and 912 MB (FastANI); the sketch store is ~5× smaller than skani's (Supplementary Fig. S6). At database scale, all-vs-all `triangle` on 5,000 GTDB-R207 representative genomes completes in ~14 min on a 32-core node; `search` against the same database takes ~4 min for 1,000 queries. skani remains faster at scale, a trade-off we state openly.

---

## Discussion

Syn2bANI demonstrates that fixed restriction-site anchors can support accurate, interpretable ANI estimation while additionally yielding synteny and structural information in the same computational pass. The chain-restricted likelihood cleanly separates divergence from shared content, and the calibration layer corrects the rate-heterogeneity bias that simulations alone fail to predict. Two boundaries are worth stating. First, skani remains faster for very large all-vs-all screens and is more accurate in the near-clonal 97–100% regime; the hybrid rule lets Syn2bANI recover the best raw estimate there. Second, the v5 calibration is trained on complete genomes and should not be silently applied to MAGs — the raw gated estimate is preferred for draft inputs, and the consistency flag marks pairs at the identifiability boundary.

The structural outputs are not an afterthought: they are a direct consequence of using coordinate-carrying anchors. In an era of strain-resolved metagenomics, ANI alone is insufficient to decide whether two genomes are truly equivalent. The synteny score and SV calls provide a second axis of comparison that is orthogonal to divergence and biologically interpretable.

---

## Online Methods

### In-silico digestion and IUPAC-aware tag geometry

Each genome is digested in silico with every enzyme in the panel (default BcgI, AlfI, AloI, FalI; all 16 Type IIB enzymes of the 2bRAD-M panel are supported [8,19]). Recognition patterns are decomposed into fixed anchor bases and IUPAC degenerate positions; matching uses precomputed 4-bit base masks so each degenerate check is one load and one bitwise AND. The digestion enforces the *complete* biological site on both strands — during validation we found that the patterns for HaeIV, Hin4I, and BaeI had silently dropped the trailing fixed anchor bases of the degenerate-containing right anchor (~16× excess tag density for HaeIV), and this was fixed alongside the likelihood geometry. Tags longer than 32 bp (CspCI, 33 bp) exceed the 2-bit packing width and are refused with a warning rather than silently mismatched — a policy instituted after a truncation asymmetry was found to destroy 91% of reverse-strand anchors for such enzymes.

For the likelihood, each enzyme carries a site geometry `{tag_len, exact_site, d2, d3}` parsed from its IUPAC anchor string: `exact_site` fully specified positions (a mutation deletes the tag), `d2` two-of-four degenerate positions (survive a mutation with probability a + (1−a)/3), `d3` three-of-four positions. The homogeneous likelihood uses the exact convolution over degenerate classes; the heterogeneous likelihood uses a negative binomial with effective body b + d2/3 + 2·d3/3 plus a closed-form survival factor. With the fully specific default panel (d2 = d3 = 0) both paths are bit-identical to the pre-fix code. Tags are stored strand-canonically (the lexicographically smaller of a tag and its reverse complement), so contig orientation is immaterial; reverse-complement self-controls on real drafts return 99.9999.

### Anchoring, chaining, and the adaptive chain-break test

Reference tags are indexed by packed sequence; sequences occurring more than `max_occurrence` (default 5) times in either genome are dropped as repeats/paralogs. Tolerant seeding (mismatch budget `tol` = 2 by default) uses the pigeonhole principle: each tag is split into `tol + 1` parts and any pair within `tol` mismatches shares at least one part exactly, so seed recovery needs no whole-genome scan. The 2-mismatch budget is what keeps the method alive below ~93% ANI — at 90% ANI a 32 bp tag matches exactly only 3.4% of the time.

Anchors are grouped by (query contig, reference contig, orientation) — never across contigs, since an inter-contig "gap" joins unlinked sequences — and chained by gap-penalized collinear DP with a bounded predecessor window. The chain-break test counts *skipped query tag positions*, not base pairs, because a fixed bp threshold provably cannot work across the ANI range: at 95% ANI a tag anchors with probability ~0.79 (nine consecutive failures are already implausible, ~6 kb), while at 85% the probability is ~0.12 (a hundred consecutive failures are ordinary, ~78 kb). Measured on exact-truth sweeps, a permissive 50 kb bp-limit lets chains bridge accessory blocks (ANI biased down 0.74 points, AF inflated), while a strict 10 kb limit fragments chains and drops poorly anchored regions (ANI biased up to +3.07 at 85% ANI). Counting skipped tag positions is scale-free and separates the two events a bp limit conflates: a deletion removes reference tags while surviving query tags stay adjacent (zero skipped positions, no break), whereas a length-preserving non-homologous replacement skips every query position inside it (break). The threshold is j* = ln(α)/ln(1−p) with α = 1e−6, floored at 5 because repeat-masked tags and panel coverage gaps occupy tag positions without being anchors. Pass 1 chains permissively and fits a provisional divergence; pass 2 re-chains at j* and re-fits. Inside each chain, query→reference coordinates are interpolated from the chain's own anchors and tolerant matching is attempted only within `local_window` bp of the predicted position — seed-and-extend, which cannot invent matches between positionally unrelated tags. Chains are processed largest-first and each query tag is counted at most once.

### Chain-restricted stratified MLE, gating, and flags

For a query tag inside a chained region, with tag length k, recognition-site length s, mutable body b = k − s, and mismatch budget tol, the outcome probabilities under per-base identity a are: found with m body mismatches, P_m(a) = C(b, m)·(1−a)^m·a^(k−m) for m ≤ tol; not found, P_miss(a) = 1 − Σ P_m(a). Mismatches are observable only in the body — a site mutation deletes the tag — so the site contributes a^s but never appears in the histogram; this also sets a hard retention ceiling (~0.72 at 95% ANI, ~0.50 at 90%). The log-likelihood is summed over per-enzyme strata (each with its own k and b) and maximized over the scalar a by golden-section search. Two partial estimators — `ani_from_loss` (miss rate only) and `ani_from_hist` (histogram only, renormalized by P(found)) — are independent estimates of the same quantity and power the diagnostics.

The heterogeneous model gives each region a rate multiplier r ~ Gamma(α, α); because rate variation acts at kilobase scale while a tag is ~30 bp, all sites in a tag share one r, and integrating r out turns the mismatch count negative binomial with two parameters (mean divergence d, shape α) identifiable from the three histogram degrees of freedom plus the miss count. Guardrails: the second parameter is spent only when the likelihood-ratio statistic exceeds 3.841, and α is clamped to [0.1, 200]. Because the gamma shape and mean couple at the identifiability boundary when few tags survive (overshoots of 4–10 ANI points at mid-ANI), the shipped point estimate is **gated**: report the homogeneous fit when |ani_from_loss − ani_from_hist| > 5 ANI points, else the gamma fit. The threshold is an effect size, not a significance level — significance-scaled variants invert their error ranking on GTDB — and was chosen at the flat optimum (4.5–6 points) of a threshold sweep on the GTDB-ANIm matrix. It fires on 5.5% of GTDB pairs (MAE 2.819 vs 2.881 always-gamma; per-pair oracle bound 2.788), 12/15 mid-ANI pairs (4.482 → 0.959), and never on uniform-rate sims, mosaic sims, or oral/gut same-species pairs, preserving gamma's advantage exactly where it is real. The `flag` column reports `BELOW_DETECTION` (expected retention < 0.20; precedence), `INCONSISTENT` (the gate fell back, or the chains carry > 0.5 rearrangement breakpoints per anchor — a structural statistic that does not share the chain-restricted likelihood denominator), else `ok`. The recalibrated flag never inverts its error ranking on any validation set (kept vs flagged MAE: GTDB-ANIm 2.415 vs 4.153; oral/gut 0.514 vs 4.269; mid-ANI 0.343 vs 1.113), at a disclosed cost in near-clonal sensitivity.

### Short-contig rescue pass

A query contig with < 8 × `min_chain_anchors` in-panel tags (~25 kb with the default panel) and no accepted chain is folded into the likelihood without chaining if it carries a collinear group of ≥ 2 unique anchors (unique = occurring once in each genome; at least one anchor with ≤ 1 mismatch, since a solely-2-mismatch basis can be a paralog; ties break on total mismatch count). Counting replicates exactly what pass-2 chaining would count on the intact genome: tags between the outer basis anchors are counted (hits to the histogram, no match to the miss count), and a tail tag beyond an outer anchor is counted only while a bracketing anchor from any contig follows within (`max_skip` + 1) × mean tag spacing along the reference. Rescued contigs contribute AF spans with the same half-median-gap extension rule as chains but create no chains, so synteny and SV outputs are untouched. The pass is gated on contig tag count, so complete and high-N50 assemblies take a bit-identical path. The residual +0.20 drift at N50 5 kb comes from contigs with 0–1 placeable anchors (~1.3 Mb on the Sakai ladder), which carry no positional evidence; rescuing them would require guessing homology.

### Database-scale screen

Batch subcommands (`dist`, `search`, `triangle`, `db search`) run a two-stage pipeline. Stage 1 reduces each tag to one strand-canonical key — its centered 18 bp packed window — per enzyme; a pair passes iff shared keys ≥ 3 AND shared/min(smaller key set) ≥ 0.001. At 18 bp, generic background homology (shared genes at 65–75% identity) no longer saturates key containment, and the gate was calibrated on GTDB-R207 to false-reject 0/500 validated ≥ 80% ANI pairs while rejecting ~83% of random pairs (measured false-reject rate on estimator-reportable pairs at n = 2,000: 1/844). In the FRR-critical 80–85% band the weakest validated true pair clears the floors by 10×/6×. Stage 2 refines survivors with the identical `chain_ani` code path and row formatter as the pairwise `ani` command, so database output is byte-identical to pairwise output on the same pairs. Screen thresholds are panel-specific and exposed as flags.

### Ridge calibration protocol

The calibration model is ridge regression on nine Syn2bANI-internal features — `ani_gated`, `ani_uniform`, `af_query`, `af_reference`, `std_err`, `retention`, `n_anchors`, `n_chains`, `n_tags_in_chains` — with median imputation and standardization, regularization α = 10 selected by inner cross-validation (scikit-learn [18]), trained on 2,520 ANIm-truth pairs (2,053 finite pairs from the stratified 2,074-pair GTDB benchmark plus 467 targeted 95.0–99.15% ANIm pairs from a skani screen of all 2.26 M same-genus GTDB-R207 pairs followed by dnadiff on 650 selected candidates; 607 species, 26 phyla, zero overlap with the existing set). Validation is **band-holdout**: each of the four ANI bands (80–85, 85–90, 90–95, 95–99) is predicted by a model trained on the other three only. Training rows are shuffled with a fixed seed before CV (the inner KFold of the ridge CV is order-sensitive); seed spread over 5 seeds is ≤ 0.0015 MAE. An expanded 18-feature variant and a gradient-boosted model on the same features were evaluated and rejected (nonlinearity buys nothing at n ≈ 2,500, and the expanded set degrades out-of-distribution on near-clonal pairs). A learning-curve analysis (492/984/1,476/1,969 training pairs → CV MAE 0.952/0.902/0.887/0.893) shows the linear model plateaus at ~1,500 pairs. The model ships as a versioned JSON swapped into the binary; all feature matrices were re-run with the post-rescue binary before training, so there is no version skew between the shipped estimator and the deployed model. Calibrated output is a separate column (`--calibrate`); it returns no value on `BELOW_DETECTION` pairs rather than extrapolating.

### The spatial-model negative result (identifiability analysis)

The mechanistic alternative to calibration was evaluated in full (prototype in Python against 19 exact-truth mosaic simulations, then the GTDB-ANIm matrix). The bias decomposes into an in-chain distribution-family term — asymptotic in-chain MAE 2.25 for gamma, 1.25 for a capped grid NPMLE, so a nonparametric fit fixes at most ~0.1–0.2 GTDB MAE — and a coverage term: the divergence of the unchained fraction is not identifiable from tag data, since the out-of-chain anchor residual is ≈ 0 within multi-match noise whether the unchained mass is saturated-divergent or accessory, and an oracle given the exact identity of the chained sample still errs by MAE 1.36 against whole-genome truth. No candidate (AF-weighted mixtures, discrete mixtures, ascertainment-aware tilts, capped NPMLE with and without LRT gating, model averaging) beats the gated baseline while holding all simulation gates. The raw estimator is therefore at its identifiability floor at 4-enzyme tag density; this analysis is why calibration, not a new likelihood, is the deployed correction layer.

### Structural-variant calling

The `struct` subcommand consumes the chains and anchors of a pairwise comparison. Inversions are orientation flips between consecutive chains; indels are coordinate gaps between consecutive chains exceeding `indel_min` (default 1 kb), with a guard requiring that no third chain anchors inside the spanned reference interval (without it, a small relocated block sandwiched between collinear chains is misreported as a deletion sized by the reference jump — observed as spurious 42–961 kb calls before the guard); translocations are chains violating global reference order. Resolution is limited by tag spacing and by unphaseable repeats at junctions (observed endpoint offset ≤ 5.5 kb on the 779 kb W3110 inversion); events whose flanking chains have < 4 anchors are not called, and equal-length homologous replacements that keep tags phased are invisible by design. Validation against dnadiff used span-based one-to-many matching because dnadiff fragments large accessory regions into many small events while Syn2bANI calls each junction once.

### Simulation framework with exact truth

All exact-truth simulators evolve sequence from a public reference (*E. coli* K-12 MG1655, ENA U00096.3, 4,641,652 bp) by applying a counted number of substitutions, so true ANI = 1 − n_subs/L exactly. Families: an ANI ladder (12 levels, 85–99.9%, each with a 400 kb inversion, with and without ~46 deletions of 200–2,000 bp); an indel sweep (true ANI 95.000; 0–4 deletions/100 kb); fragmentation (20–201 contigs, no sequence lost, ~50% reverse-complemented, order shuffled); accessory confounds (core ANI 95.000; 0–50% of the genome replaced by composition-preserving shuffled blocks, plus a block-count control at fixed total fraction); mosaic/rate-heterogeneity cases (per-block rates ~ Gamma(α, α), α = 0.5/1/2, mean ANI 90–98, plus deliberately misspecified bimodal cases); and GC coverage (ladders on five real template genomes, GC 27.2–72.1%). Deletions do not change true ANI (counted substitutions over reference length), so the deleted fraction doubles as AF ground truth. Substitutions are uniform over the three alternative bases with no GC-biased mutation model, transition/transversion ratio, codon structure, or selection — a deliberate simplification justified by the finding that GC does not explain the real-data degradation and that site turnover emerges from uniform substitution on a real genome without any special mechanism. Enzyme digestion is never simulated; the released binary extracts tags itself.

### Benchmark datasets, truth, and tools

All datasets, their ground truth, and the analyses they support are summarized in Table 3. The simulation-generation framework is described in Methods 4.9.

**Sampling frame.** GTDB release R207 [11,12] is the primary sampling frame for the real-genome benchmarks; all GTDB-R207 comparisons use dnadiff/ANIm (MUMmer 3.23, 1-to-1) as truth. CAMI2 (strain-madness and marine short-read campaigns, with source genomes and per-read ground truth) is the sampling frame for the MAG benchmark.

**Simulated datasets (exact truth by construction).** Seven families, all evolved from *E. coli* K-12 MG1655 except the GC ladders: (i) the **ANI ladder** — 12 levels, 85–99.9%, each genome carrying a 400 kb inversion, generated with and without ~46 deletions of 200–2,000 bp; (ii) the **indel sweep** — ANI 95.000, 0–4 deletions per 100 kb; (iii) the **fragmentation ladder** — 20–201 contigs, ~50% reverse-complemented, order shuffled, no sequence lost; (iv) the **accessory-fraction sweep** — core ANI 95.000 with 0–50% of the genome replaced by composition-preserving shuffled blocks, plus a block-count control at fixed total fraction; (v) the **mosaic/rate-heterogeneity family** — per-block rates ~ Gamma(α, α), α = 0.5/1/2, mean ANI 90–98, plus deliberately misspecified bimodal cases (19 used in the identifiability analysis of Methods 4.7); (vi) the **GC ladders** — substitution ladders on five real template genomes spanning GC 27.2–72.1%; and (vii) the **inversion ladder** for synteny validation — ANI 95.00/98.00 with 0–32 non-overlapping inversions of 100–400 kb (2 exact breakpoints each; `scripts/synteny_bench/`).

**GTDB-R207 benchmark series.** These datasets form a single coherent evaluation built from GTDB-R207. They are disjoint at the genome level wherever independence is required.
- **Calibration/training set (v5).** 2,520 pairs (80–99.5% ANIm) used to train the ridge calibration. This combines 2,074 representative pairs stratified across 80–85/85–90/90–95/95–99 bands and phyla, with 467 targeted 95–99.5% pairs selected by screening all 2.26 M same-genus representative pairs with skani and running dnadiff on 650 candidates (607 species, 26 phyla).
- **Large-scale FastANI comparison set.** 45,000 representative pairs sampled across taxonomic levels (intra-species / intra-genus / intra-family / random; `scripts/sample_gtdb_pairs.py`, seed-fixed). FastANI reports 672 of these pairs and they provide the large-scale comparison in Fig. 6.
- **Strict held-out benchmark.** 43,334 same-genus representative pairs sampled from the 628,617 non-self GTDB-R207 representative pairs passing a skani pre-screen (AF ≥ 15%), stratified into bands 80–85 / 85–90 / 90–95 / 95–100% (12,172 / 16,000 / 14,758 / 404 pairs; 24,831 genomes, 74 phyla; `scripts/gtdb50k/`). Every genome in the v5 calibration set was excluded in both directions, making this a strict holdout for the gate, flag, and calibration; dnadiff/ANIm truth was computed for all pairs.
- **High-ANI test set.** 727 test pairs sampled from non-representative GTDB-R207 genomes (95–97%, n = 95; 97–100%, n = 632), excluding genomes in the v5 calibration set or the 43,334 held-out set in either direction. These add dense coverage of the 95–100% regime that the held-out representative set lacks.
- **Unified 80–100% benchmark.** The 43,334 held-out pairs combined with the 727 high-ANI test pairs (44,061 total), used to evaluate the hybrid estimator (Supplementary Fig. S3).

**External validation sets.** Independent of the GTDB-R207 series.
- **Mid-ANI set.** 15 pairs (*Bifidobacterium longum* × *B. breve*, *Veillonella parvula* × *V. dispar*), ANIm truth 87.6–90.2%; ANIm alignment coverage is only 46–72% of these genomes.
- **Oral/gut set.** 50 isolate genomes from 10 species — 5 oral (*Aggregatibacter actinomycetemcomitans*, *Fusobacterium nucleatum*, *Porphyromonas gingivalis*, *Streptococcus mutans*, *S. sanguinis*) and 5 gut (*Akkermansia muciniphila*, *Bacteroides fragilis*, *B. longum*, *Faecalibacterium prausnitzii*, *Ruminococcus intestinalis*), 5 genomes per species — compared all-vs-all (1,225 pairs). The 100 same-species pairs serve as the external near-clonal benchmark against FastANI and skani; 1,100 cross-species pairs probe the below-detection regime.

**Structural-variation and MAG validation.**
- **Enterobacteriaceae completes.** 13 finished chromosomes versus *E. coli* K-12 MG1655, spanning ANI ~81–99.99%. Three pairs (MG1655/W3110, MG1655/O157:H7 Sakai, *Salmonella* Typhi CT18 / *S.* Typhimurium LT2) double as the structural-variation validation set with dnadiff structural truth.
- **Real draft assemblies.** 8 *E. coli* assemblies from ENA (GCA_001283865, GCA_001077875, GCA_001284645, GCA_001283245, GCA_001283605, GCA_001284145, GCA_001283205, GCA_001075925; 88–8,025 contigs), including reverse-complement self-controls.
- **CAMI2 MAG benchmark.** 35 short-read samples (25 strain-madness, 10 marine) assembled with MEGAHIT, binned with MetaBAT2, and quality-scored with CheckM2; 695 bins ≥ 100 kbp were retained (HQ 200 / MQ 210 / LQ 285). Contigs were assigned to CAMI2 source genomes by exact k-mer membership; each bin was classified clean, strain-mixed, or cross-species, and paired with its dominant source genome ("anchor"). Each tool estimated ANI of the bin against the anchor reference; dnadiff/ANIm provided truth. Pipeline scripts and per-pair outputs are in `results/mag_validation/`.
- **Syntracker validation isolates.** Four published isolate collections were used to test whether high ANI guarantees conserved genome architecture: *Escherichia coli* hypermutator (23 isolates, 253 pairs), *Helicobacter pylori* (77 isolates, 2,926 pairs), *Neisseria gonorrhoeae* (66 pairs), and *Streptomyces rimosus* (190 pairs). SRA accessions and references are in `data/syntracker/`; all isolates were compared all-vs-all with `syn2bani ani`, and `synteny_score` and `breakpoint_count` were taken from the same pass. Top-discordant cases were re-analyzed with `syn2bani struct` (Fig. 11, Fig. 12; `results/syntracker_validation/`).

**Tools and versions.** Syn2bANI 0.1.0 (default panel BcgI,AlfI,AloI,FalI; deployed calibration v5), skani 0.1.0 (0.3.2 for the HPC database-scale runs), FastANI 1.33/1.34, MUMmer/dnadiff 3.23. Per-pair efficiency was measured on a Mac Studio (Apple M4 Max, 16 cores, 128 GB RAM; 3 repetitions, medians reported); database-scale runs used a 32-thread HPC node. skani's n = 22 `dist` time covers only 302/484 reported pairs (minimum-AF filter) and excludes sketching; we note this wherever the comparison appears.

### Implementation

Syn2bANI is implemented in Rust (edition 2021) with needletail [13] for FASTA parsing and rayon [14] for data parallelism. The CLI mirrors skani's subcommand structure (`dist`, `search`, `sketch`, `triangle`) and adds `ani` (pairwise MLE), `db` (sketch-database management), and `struct` (SV calling). The `.s2ba` sketch format stores 2-bit-packed tags with their coordinates and records the enzyme table (older sketches remain readable); a sketch is ~120–180 KB per genome and the store is ~5× smaller than skani's at n = 22. Sketch reuse makes repeated comparisons against the same genomes bit-identical at ~5× lower peak memory. The calibration model is an embedded JSON (ridge weights, scaling, feature order) evaluated at runtime in < 1 µs; the 95–99.5% training expansion and all calibration scripts are released. The test suite (101 tests at the time of writing, including rescue-pass, geometry, screen, and synteny-statistics regressions) is green in release mode.

---

## Data availability

All benchmark datasets, ground-truth files, evaluation matrices, and figure-generation scripts are available in the project repository at https://github.com/HuangShiLab/Syn2bANI-paper. These include the simulation manifests, the 2,074-pair GTDB-R207 ANIm benchmark (`results/panel_by_band/`), the 467-pair 95–99.5% ANIm expansion (`results/anim_truth_hi95.tsv`), the 43,334-pair held-out GTDB-R207 benchmark (`results/gtdb50k/`), the unified high-ANI benchmark (`results/gtdb50k/high_ani_results.tsv`), the CAMI2 MAG benchmark (`results/mag_validation/`), the SV validation (`results/sv_validation/`), the synteny benchmark (`results/synteny_bench/`), and the Syntracker isolate re-analysis (`results/syntracker_validation/`). GTDB release R207 is available from https://gtdb.ecogenomic.org. The reference for exact-truth simulations is *E. coli* K-12 MG1655, ENA accession U00096.3.

## Code availability

Syn2bANI is free and open source (MIT License) at https://github.com/HuangShiLab/Syn2bANI, including the Rust source, integration tests, the exact-truth simulation harness, and the embedded calibration model. The version benchmarked here is 0.1.0 with calibration model v5 (commit `fe0f36c`).

## Acknowledgements

[To be determined.]

## Author contributions

[To be determined — placeholder.]

## Competing interests

[To be determined — placeholder: The authors declare no competing interests.]

## References

1. Konstantinidis, K. T. & Tiedje, J. M. Genomic insights that advance the species definition for prokaryotes. *Proc. Natl. Acad. Sci. USA* **102**, 2567–2572 (2005).
2. Jain, C., Rodriguez-R, L. M., Phillippy, A. M., Konstantinidis, K. T. & Aluru, S. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. *Nat. Commun.* **9**, 5114 (2018).
3. Bowers, R. M. et al. Minimum information about a metagenome-assembled genome (MIMAG) of bacteria and archaea. *Nat. Biotechnol.* **35**, 725–731 (2017).
4. Marçais, G. et al. MUMmer4: A fast and versatile genome alignment system. *PLoS Comput. Biol.* **14**, e1005944 (2018).
5. Richter, M. & Rosselló-Móra, R. Shifting the genomic gold standard for the prokaryotic species definition. *Proc. Natl. Acad. Sci. USA* **106**, 19126–19131 (2009).
6. Shaw, J. & Yu, Y. W. Fast and robust metagenomic sequence comparison through sparse chaining with skani. *Nat. Methods* **20**, 1661–1665 (2023).
7. Wang, S., Meyer, E., McKay, J. K. & Matz, M. V. 2b-RAD: a simple and flexible method for genome-wide genotyping. *Nat. Methods* **9**, 808–810 (2012).
8. Sun, Z., Huang, S., Zhu, P. et al. Species-resolved sequencing of low-biomass or degraded microbiomes using 2bRAD-M. *Genome Biol.* **23**, 36 (2022).
9. Yang, Z. Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites: approximate methods. *J. Mol. Evol.* **39**, 306–314 (1994).
10. Enav, H., Paz, I. & Ley, R. E. Strain tracking in complex microbiomes using synteny analysis reveals per-species modes of evolution. *Nat. Biotechnol.* (2024).
11. Parks, D. H. et al. A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life. *Nat. Biotechnol.* **36**, 996–1004 (2018).
12. Parks, D. H. et al. GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. *Nucleic Acids Res.* **50**, D785–D794 (2022).

---

*Manuscript draft — Nature Methods format. Main-text target: ~3,000–5,000 words; abstract: 153 words.*
*Corresponding author: Shi Huang (huangshi@njau.edu.cn)*

---

## Figure Legends

**Figure 1 | The Syn2bANI estimator.** (a) In-silico digestion with the default Type IIB enzyme panel (BcgI, AlfI, AloI, FalI); tags carry sequence and coordinate. (b) Shared tags become anchors; lost tags contribute to the miss count. (c) Anchors are chained per (query contig, reference contig, orientation) by gap-penalized collinear dynamic programming in two passes; inverted blocks are detected as orientation flips. (d) Per-enzyme tag outcomes inside chains enter a chain-restricted maximum-likelihood fit; outputs include gated ANI, standard error, aligned fractions, synteny score, SV calls, reliability flag, and `ani_upper95`. Source: `figures/report/fig1_algorithm_schematic.png`.

**Figure 2 | Accuracy and robustness under exact truth.** (a,b) Twelve *E. coli* genomes evolved at true ANI 85.0–99.9% (400 kb inversion + deletions): estimated vs true ANI and signed error. MAE: Syn2bANI 0.073%, skani 0.377%, FastANI 0.742%. (c) Indel sweep at 95% ANI (0–4 deletions/100 kb): MAE 0.081. (d) Simulated fragmentation (20–201 contigs): MAE 0.093. (e) GC ladders on five genomes (27.2–72.1% GC): MAE 0.074–0.356. (f) Accessory-content sweep (0–50% shuffled blocks at 95% core ANI): ANI estimate flat, `af_query` tracks true shared fraction. Source panels: `figures/report/fig1_simulation_ladder.png` and `figures/report/fig2_robustness.png`.

**Figure 3 | GTDB-R207 benchmark against ANIm truth.** (a) Band-holdout cross-validation on 2,520 training pairs: MAE by ANI band for raw gamma, calibrated v5, skani, and FastANI. (b) Strictly held-out 43,334 same-genus pairs with training genomes excluded: calibrated Syn2bANI v5 vs skani and FastANI. Source: `figures/report/fig7_anim_by_band.png` and `figures/report/fig_gtdb50k_heldout.png`.

**Figure 4 | Unified 80–100% GTDB-R207 benchmark and hybrid estimator.** 43,334 held-out pairs plus 727 high-ANI test pairs. (a) Estimated vs ANIm truth for raw gated, v6-calibrated, hybrid (raw ≥98% / calibrated <98%), skani, and FastANI. (b) Per-band MAE. (c) Signed-error distributions. (d) High-ANI zoom. Source: `figures/gtdb_r207_unified_benchmark.png`.

**Figure 5 | Accuracy on 695 CAMI2 MAGs.** (a) Syn2bANI raw gated estimate vs dnadiff ANIm truth, colored by contamination class. (b) Absolute-error distributions by tool. (c) Syn2bANI error by CheckM2 quality tier. Source: `figures/report/mag_validation.png`.

**Figure 6 | High-ANI pairs hide extensive rearrangements.** (a) *E. coli* hypermutator and *H. pylori* isolates: ANI vs synteny score (Syn2bANI left, skani right); highlighted pairs were re-analyzed with `syn2bani struct`. (b) *N. gonorrhoeae* and *S. rimosus* isolates show the same ANI–synteny decoupling. Source: `figures/syntracker_validation/syntracker_high_ani_low_synteny.png` and `figures/syntracker_validation/syntracker_supp_ngonorrhoeae_srimosus.png`.