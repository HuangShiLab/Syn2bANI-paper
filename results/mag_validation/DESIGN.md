# Task 1 Design: Real-MAG validation benchmark for syn2bANI

Status: DESIGN (scoping done 2026-08-19, read-only). Nothing executed yet.
Goal: validate syn2bANI on **real MAGs with realistic binning error and contamination**
(>=100 genomes), against skani and FastANI, with a defensible ANIm truth.

---

## 1. What is actually on the HPC (verified 2026-08-19)

### CAMI2 datasets under /lustre1/g/aos_shihuang/data/

| Dataset | Reads | Source genomes | Gold-standard truth present? |
|---|---|---|---|
| `cami2_strain` (strain madness) | 100 samples extracted, each `short_read/sample_N/reads/anonymous_reads.fq.gz` (~2 GB gz, anonymized headers `S0R0/1`) | 408 FASTA in `source_genomes/short_read/source_genomes/` (+ .fai; includes MIKI-NS* unpublished strains) | **No** per-read mapping / GS binning downloaded |
| `cami2_marine` | 10 samples extracted (`marmgCAMI2_sample_N_reads/2018.08.15_09.49.32_sample_N/reads/anonymous_reads.fq`), ~5.2 GB tar each | `genomes/`: 977 files = 953 species-named references + 24 assembled `N_SN_contigs.fasta` (SPAdes-style headers, the unpublished strains) | **Yes**: per-sample `reads_mapping.tsv.gz` (read -> genome_id `OtuNNNN.0`, tax_id); `ground_truth_from_unidentified_source/` (species abundance matrix) |
| `cami2_plant_associated` | 21 sample tars NOT extracted; 43 tarballs | 894 FASTA in `source_genomes/` | No |

Other relevant data:
- `gtdb-r207/genomes_all/`: 65,703 GTDB r207 genomes (complete+draft reps) — usable as the "reference database" side of MAG-vs-DB queries.
- `cami2_marine/RefSeq_genomic_20190108/`: full RefSeq snapshot + ncbi_taxonomy (large; only if needed).
- `cami2_marine/GTDBtk_r207_results/`: GTDB-Tk run on the 977 marine source genomes (identify step present).

**Key finding: no ready-made MAGs/bins exist anywhere.** CAMI2 ships reads + source genomes; the MAGs
must be produced by our own assembly+binning run. This is standard (AMBER/CAMI2-style evaluation)
and is in fact a *feature*: we control the binner and can derive per-contig truth exactly.

### Tools verified on HPC (all present)

| Tool | Path / env | Version |
|---|---|---|
| syn2bani | `/lustre1/g/aos_shihuang/Syn2bANI-hi95/target/release/syn2bani` | 0.1.0 (068119c+, current) |
| skani | `/group/aos_shihuang/conda/envs/gtdbtk310/bin/skani` | 0.3.1 |
| fastANI | `/group/aos_shihuang/conda/envs/fastani/bin/fastANI` (also in anvio, coverm, gtdbtk*) | present |
| dnadiff / nucmer | `/group/aos_shihuang/conda/envs/anvio/bin/` | dnadiff 1.3 (MUMmer) |
| MEGAHIT | `/group/aos_shihuang/conda/envs/megahit/bin/megahit` | v1.2.9 |
| MetaBAT2 / MaxBin2 / CONCOCT | conda envs `metabat2`, `maxbin2`, `concoct` | present |
| CheckM2 (+DB) / CheckM1 | env `checkm2`; DB at `/lustre1/g/aos_shihuang/db/CheckM2_db` | 1.1.0 |
| coverm (with minimap2), bowtie2, seqkit | envs `coverm`, `bowtie2`, `seqkit` | present |

Slurm: `intel` (84 nodes) and `amd` (63 nodes) partitions, 7-day limit. Nothing in user PATH;
always call binaries by absolute path.

---

## 2. Recommended dataset(s)

**Primary: CAMI2 strain madness** (`cami2_strain`). Rationale: 408 source genomes with heavy
strain-level redundancy (the hardest and most paper-relevant case for an ANI estimator);
100 samples give ample MAG yield from a subset.

**Secondary: CAMI2 marine** (`cami2_marine`). 10 samples, 977 source genomes, different biome,
and it is the only dataset with the per-read truth mapping already on disk (useful for
cross-checking the contig-assignment pipeline).

Plant-associated: skip (tars not extracted; adds nothing over the other two).

---

## 3. Cohort definition

1. Assemble + bin a subset: **25 strain-madness samples** (random, seed 42) + **all 10 marine samples**.
2. Bin with **MetaBAT2** (default; optionally MaxBin2 on 5 samples as a binner-robustness check).
3. Keep bins >= 100 kbp total. Expected yield: ~30-80 bins/strain sample, ~50-150/marine sample
   => well over 1,000 raw bins; plenty for a 100+ MAG cohort.
4. Per-contig truth assignment: map every bin's contigs to the dataset's source genome set
   (minimap2 `asm5` from `coverm` env, or nucmer from `anvio`). Assign a contig to a source genome
   if best hit covers >= 80% of contig length at >= 95% identity; else "unassigned".
   Cross-check marine assignments against `reads_mapping.tsv.gz`-derived expectations.
5. Derived per-bin quantities:
   - `majority_genome`: source genome with most assigned bp.
   - `completeness_est` = covered fraction of majority genome (alignment-based), plus CheckM2
     completeness/contamination (DB present) for MIMAG tiers:
     HQ (>=90% comp, <=5% cont), MQ (>=50%, <=10%), LQ (rest).
   - **Contamination classes from gold standard** (alignment-based, independent of CheckM2):
     (a) clean: >=95% assigned bp from majority genome;
     (b) strain-mixed: contaminant bp from other strain(s) of the *same species*;
     (c) cross-species contaminated: contaminant bp from other species. Report contamination % in bp.
   - N50, #contigs, total size (fragmentation metrics, tie into the fragmentation-ladder narrative).
6. Final cohort: all bins, but analyses reported for (i) HQ+MQ clean, (ii) strain-mixed,
   (iii) cross-species contaminated, (iv) LQ. Target >=100 MAGs in the HQ/MQ core cohort —
   comfortably met.

---

## 4. Truth strategy (ANIm)

For each MAG, compute **dnadiff** (MUMmer, `/group/aos_shihuang/conda/envs/anvio/bin/dnadiff`)
against its `majority_genome` source FASTA. Record dnadiff's 1-to-1 ANI **and** aligned bases
(AF) for both directions. Accept a pair as truth if AF >= 60% on the MAG side; for AF 30-60%
keep but flag as low-AF; below 30% report separately (this is where all tools should degrade —
it is part of the result, not noise).

Truth subset size: dnadiff ~8-10 s/pair. Run **all** MAG–majority-genome pairs (expected
~600-1,200 pairs => ~2-4 core-hours, trivial in an array job), so no sampling bias; if yield
explodes, subsample ~500 stratified by quality tier x contamination class.

Second truth axis (reference-DB scenario): each MAG vs its **nearest GTDB r207 representative**
(from `gtdb-r207/genomes_all`, picked by skani pre-screen). No independent truth exists for
MAG-vs-rep beyond the tools themselves, so this axis reports *cross-tool agreement* and flag
behavior rather than accuracy vs truth. For strain-madness MAGs whose source genome is itself
a real (RefSeq) genome, the MAG-vs-source dnadiff truth still anchors absolute accuracy.

Sanity: FastANI vs dnadiff concordance on the truth subset (expected r>0.99 from prior GTDB work)
to catch truth-computation bugs before interpreting syn2bANI.

---

## 5. Evaluation metrics

- MAE, median bias, RMSE, within-{0.1,0.5,1.0}% ANI rates: per tool x quality tier x contamination class.
- Error vs N50 / #contigs (degradation curve; compare with simulation ladder).
- Error vs alignment-based contamination % (the key new axis).
- Contamination-flag behavior: does syn2bANI's gating/contamination flag fire on class (b)/(c)
  bins? Report flag precision/recall vs gold-standard contamination at 1%, 5%, 10% bp thresholds,
  and whether `ani_upper95` brackets the truth.
- AF-aware analysis: error vs dnadiff AF; confirm errors concentrate in low-AF pairs.
- Cross-tool table: syn2bANI vs skani 0.3.1 vs FastANI on identical pairs (all three are
  sketch/minimap-fast: run on the full cohort, not just the truth subset).

---

## 6. Execution plan (slurm, HPC workdir `/lustre1/g/aos_shihuang/mag_validation/`)

| Job | What | Resources | Est. wall |
|---|---|---|---|
| J1 assembly | MEGAHIT per sample, array 35 (25 strain + 10 marine) | 16 cpu, 64G, intel | ~6-12 h |
| J2 depth | bowtie2 map reads->contigs + coverm/jgi depth, array 35 | 16 cpu, 32G | ~3-6 h |
| J3 binning | MetaBAT2 per sample (>=2.5kb contigs), array 35 | 8 cpu | <1 h |
| J4 quality | CheckM2 predict on all bins (DB present) | 16 cpu | ~1-2 h |
| J5 assign | minimap2/nucmer contigs->source genomes; build truth tables | 16 cpu | ~2 h |
| J6 fast tools | syn2bani sketch+dist / skani dist / fastANI: MAG vs majority source, top-3 same-species alternates, nearest GTDB rep | 16 cpu | <1 h (ms/pair) |
| J7 truth | dnadiff MAG vs majority source, array (~600-1200 pairs, 20/task) | 8 cpu | ~1-2 h |
| J8 collect | cat-over-ssh TSVs back to Mac (scp disabled) | login node | minutes |

Total elapsed: ~2-3 days including queueing. All outputs TSV; analysis/plots pulled back to
`results/mag_validation/` in this repo:
```
results/mag_validation/
  DESIGN.md                 (this file)
  bins.tsv                  (bin, sample, dataset, size, N50, checkm2_comp/cont, gs_comp, gs_cont, class)
  truth_dnadiff.tsv         (bin, majority_genome, anim, af_mag, af_ref)
  ani_fast_tools.tsv        (bin, reference, syn2bani, skani, fastani, flags, ani_upper95)
  figures/                  (tier-wise scatter, error-vs-N50, error-vs-contamination, flag ROC)
  MAG_VALIDATION_REPORT.md
```

---

## 7. Risks / unknowns

- **No per-read mapping for strain madness on disk.** Not blocking: contig->genome truth comes
  from alignment to the 408 source genomes (standard practice). The CAMI2 official GS binning
  applies to *their* gold-standard assembly, not ours, so alignment-derived truth is needed anyway.
  Optionally the small `reads_mapping`/GSA files can be fetched from frl.publisso.de for a
  cross-check (a few hundred MB max — needs user approval; not downloaded during scoping).
- **Strain-level truth ambiguity**: with multiple near-identical strains, a contig may hit
  several genomes at >=95% id. Mitigation: report "assignable fraction"; treat same-species
  multi-hits as that species (strain-mixed class already handles this).
- **dnadiff on LQ bins** can give very low AF — handled via AF tiers in Sec. 4 rather than discarding.
- **Binner dependence**: results could be MetaBAT2-specific; mitigated by the 5-sample MaxBin2
  check and by marine cross-validation.
- **Assembly scale**: 25+10 samples chosen to stay within ~1 day of compute; all 100 strain
  samples is a trivial array scale-up if reviewers want more.
- **Fallback if CAMI2 proves unusable**: (a) simulate MAG-like degradation from GTDB r207 draft
  genomes + synthetic contamination mixes from `genomes_all` (fully local, no download);
  (b) UHGG/MGnify real MAGs — but large downloads and no per-contig truth, so strictly worse
  than CAMI2 for contamination ground truth; (c) extend the existing ENA draft E. coli set
  (already used) — too narrow alone.
- **CheckM2 DB** confirmed at `/lustre1/g/aos_shihuang/db/CheckM2_db`; set
  `CHECKM2DB` env var accordingly in J4.

## 8. Open questions for the user

1. Sample budget: is 25 strain + 10 marine samples OK, or go straight to all 100 strain samples?
2. Download the small CAMI2 strain GS/mapping files from publisso as a cross-check (needs ~100s of MB)?
3. Include the MAG-vs-GTDB-rep axis (Sec. 4, cross-tool agreement only) or keep strictly truth-anchored?
4. Binner: MetaBAT2 only, or MetaBAT2 + MaxBin2 (+CONCOCT) for a CAMI-style ensemble (more realistic MAGs, more compute)?
