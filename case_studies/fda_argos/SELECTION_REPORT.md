# FDA-ARGOS case study — collection selection report

**Date of download:** 2026-09-25
**Source:** NCBI BioProject PRJNA231221 (FDA-ARGOS), assembly database queried via eutils (`esearch`/`esummary`/`efetch`); FASTAs via `datasets` CLI v18.33.1.
**skani:** v0.3.2 for the final gate and pair list (binary at `/Users/macstudio/.cargo/bin/skani`);
see "skani version re-check" below — the initial gate was run with v0.1.0.

## Chosen collection: *Escherichia coli* (129 genomes)

Chosen because it is the only priority candidate passing the ≥99.9-ANI decision gate (see below).
Per the task's priority order, *Salmonella enterica* was tested first and failed the gate;
*E. coli* was tested second and passed, so *Klebsiella pneumoniae* gate is reported only for completeness.

## Isolate counts per candidate species (in PRJNA231221)

Counts after deduplication by accession root (GCA/GCF mirror and old assembly versions of
the same record removed — an earlier strain-name dedup left same-assembly version duplicates
such as `GCF_001558355.1/.2` labelled "LT2"/"FDAARGOS_94", which inflated trivial 100%-ANI pairs):

| Candidate | assemblies | unique genomes downloaded | notes |
|---|---|---|---|
| *Salmonella enterica* (all serovars) | 47 | 36 | ≥15 serovars; largest same-serovar block: Typhimurium (7 strains). No Typhi, 2 Heidelberg, 0 Enteritidis. |
| *Escherichia coli* | 136 | 129 (+1 no FASTA) | largest FDA-ARGOS species block; incl. O157:H7 (4) and O104:H4 (≥3 serotyped) isolates; no ST131 labels available |
| *Klebsiella pneumoniae* | 54 | 49 | too few; 1 ST258-typed strain only per BioSample annotations |

## Decision gate (skani triangle, threads=8, measured)

| Collection | n genomes | total pairs | ANI ≥ 99.9 | ANI ≥ 99 | ANI ≥ 97 |
|---|---|---|---|---|---|
| *S. enterica* | 36 | 630 | **43 (FAIL, <50)** | 113 | 561 |
| *E. coli* | 129 | 8,256 | **102 (PASS, ≥50)** | 490 | 5,087 |
| *K. pneumoniae* | 49 | 1,176 | 11 (fail) | 664 | 1,176 |

### skani version re-check (2026-09-25)

The gate above was first run with skani **0.1.0** (the stale crates.io release that was
installed at `~/.cargo/bin/skani`; `cargo install skani --locked` only offers 0.1.1 there).
skani **0.3.2** was then built from source (`cargo install --git
https://github.com/bluenote-1577/skani --tag v0.3.2 --locked`) and the triangle re-run
(`skani_triangle_ecoli_v032.tsv`; note 0.3.2 writes a lower-triangle file, 0.1.0 wrote
upper-triangle). Result: **103 pairs ≥ 99.9 with 0.3.2 vs 102 with 0.1.0** — 101 pairs
shared, 1 old-only, 2 new-only; ANI on shared pairs shifts by ≤ 0.07 (mean 0.002).
The gate verdict is unchanged (PASS). `skani_triangle.tsv` and `pairs99.9_ecoli.tsv`
now hold the 0.3.2 results (103 pairs; the older 0.1.0 E. coli outputs are kept as
`skani_triangle_ecoli.tsv`/`.af` for provenance).

*E. coli* gate details: 103 pairs ≥ 99.9% ANI (skani 0.3.2) involving 66 of the 129 genomes;
all pairs are between distinct strain names (no technical replicates). Tightest sub-cluster:
FDAARGOS_1238/1239/1240 + FDAARGOS_1346–1356 (ANI 99.97–100.00, a QC/reference series, USA:MD)
plus a clinical O104:H4 cluster (FDAARGOS_348/400–403/772, stool/bloody-diarrhea isolates,
Georgia 2009 outbreak era). See `pairs99.9_ecoli.tsv`.

## Metadata availability

From BioSample records (129/130 fetched; full dump in `biosample_attributes_full.tsv`):

- **Available at ~full coverage:** strain, isolate, culture_collection, collected_by,
  collection_date (129/130), geo_loc_name (129/130), lat_lon, isolation_source (129/130),
  host, host_age/sex, host_disease, host_disease_outcome/stage, host_health_state.
  (Values are sometimes "missing"/"not applicable" for the QC-series isolates.)
- **Sparse:** serotype only 9-10/130 (O157:H7 and O104:H4 identifiable).
- **AMR genotype: AVAILABLE** — not in BioSample, but obtained 2026-09-25 from the NCBI
  Pathogen Detection Isolates Browser backend API
  (`POST https://www.ncbi.nlm.nih.gov/pathogens/pathogens-srv/?action=solr2txt`, query DSL
  `[display()].from(isolates).usingschema(/schema/pathogen).matching(q=="biosample_acc:(...)")`).
  127/130 accessions matched a Pathogen Detection record and **all 127 carry AMRFinderPlus
  results**: `AMR_genotypes` (e.g. blaTEM-1 36/130, sul2 35, gyrA_S83L 46), plus
  `virulence_genotypes` and `stress_genotypes` for the same 127. AST (phenotype) data:
  none for these isolates. Merged table: **`metadata_ecoli.tsv`** (raw dump:
  `pd_isolates_raw.tsv`). `epi_type` says 124/130 are clinical. The usable metadata axes
  for the paper: AMR genotype (primary), virulence genotype, epi_type, geography.

## Files

- `accessions_used.tsv` — 130 selected *E. coli* accessions with metadata; `genome_file`
  column marks FASTA presence (`MISSING` for `GCA_000783655.1`, see caveats).
- `genomes/ecoli/*.fna` — 129 FASTAs (chosen collection).
- `genomes/salmonella/`, `genomes/klebsiella/` — gate-test sets, kept on disk.
- `skani_triangle.tsv` — skani **0.3.2** triangle output for the chosen 129-genome
  *E. coli* set (= `skani_triangle_ecoli_v032.tsv`; `skani_triangle_ecoli.tsv` is the
  old 0.1.0 run kept for provenance; species-specific files and `.af` matrices retained).
- `metadata_ecoli.tsv` — merged BioSample + Pathogen Detection metadata incl. AMR genotypes.
- `struct_pairs99.9/` + `struct_summary.tsv` + `ANALYSIS.md` — syn2bani struct results
  for the 103 ANI ≥ 99.9 pairs (see ANALYSIS.md).
- `pairs99.9_{salmonella,ecoli,klebsiella}.tsv` — extracted ≥99.9-ANI pair lists
  (`ecoli` = skani 0.3.2, 103 pairs; others 0.1.0).
- `all_argos_assemblies.tsv` — full 1,707-assembly PRJNA231221 inventory.
- `tabulate_argos.py`, `fetch_ecoli_metadata.py`, `run_struct_pairs.sh` — reproducible
  fetch/analysis scripts.

## Caveats

1. FDA-ARGOS is a *reference-quality* collection, not clone-level dense sampling: the
   200–500-genome target within one species is unattainable (max 136 *E. coli* assemblies).
   The *E. coli* set meets the ANI gate but is smaller than the nominal target range.
2. `GCA_000783655.1` (strain FDA_MicroDB_82) is listed in the assembly db but its datasets
   package contains no genomic FASTA (2014 contig-level record); excluded from skani runs.
3. The tightest *E. coli* cluster (FDAARGOS_1346–1356/1238–1240) has sparse BioSample
   metadata ("missing", host "not applicable") — likely a reference/QC series rather than
   clinical isolates; the O104:H4 and O157:H7 subsets carry genuine clinical metadata.
4. skani was run with default (learned ANI) mode; ANI values are model-adjusted.
5. No FDA-hosted isolate metadata table (fda.gov "sample table") was found in a
   usable machine-readable form; BioSample records were used instead (this is the same
   metadata FDA-ARGOS submits to NCBI, 100% core-field coverage per the FDA-ARGOS paper,
   [Sichtig et al. 2019](https://www.biorxiv.org/content/10.1101/482059v1.full)).
