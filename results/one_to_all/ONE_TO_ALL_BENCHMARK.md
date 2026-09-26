# One-to-all search benchmark at full GTDB-R207 scale

Date: 2026-09-26 (final; all queries completed).

Setup: 7 case-study query genomes (GTDB R207 representative accessions:
E. coli GCF_003697165.2, B. longum GCF_000196555.1, H. pylori GCF_900478295.1,
S. aureus GCF_001027105.1, P. aeruginosa GCF_001457615.1, S. rimosus
GCF_008704655.1, N. gonorrhoeae GCF_003315235.1) searched against all 65,703
GTDB R207 representative genomes. Threads = 16 unless noted. Each run records
wall time (python timer + /usr/bin/time -v), peak RSS (time -v + /proc tree
poll), tool version, host and timestamps; see `one_to_all_benchmark.tsv`.
Runner: `scripts/bench_gtdb_one_to_all.py` (+ loop slurm scripts).

## ANI search: wall time per query (seconds)

| query | skani dist (FASTA) | skani search (sketch DB) | syn2bani dist (FASTA) | syn2bani search (sketch DB) | fastANI (FASTA) |
|---|---|---|---|---|---|
| E. coli | 1076.4 | 1308.8 (cold) | 15071.9 | 8193.0 (cold) | OOM |
| B. longum | 1105.7 | 8.8 | 14714.1 | 5886.5 | OOM |
| H. pylori | 1085.2 | 7.0 | 12983.3 | 5414.5 | OOM |
| S. aureus | 1068.5 | 43.3 | 12777.3 | 5536.1 | OOM |
| P. aeruginosa | 1030.8 | 15.6 | 13643.1 | 5536.0 | OOM |
| S. rimosus | 899.2 | 44.8 | 14749.4 | 5307.3 | OOM |
| N. gonorrhoeae | 866.3 | 6.8 | 13177.6 | 5051.1 | OOM |

## Peak RSS (GB)

| mode | peak RSS |
|---|---|
| skani dist (FASTA) | 51.6 |
| skani search | 4.3–4.5 |
| syn2bani dist (FASTA) | 51.5 |
| syn2bani search | 44.0 |
| fastANI | >130 GB at t=16, >300 GB at t=8, >511 GB at t=32, >65 GB at t=1/64 GB — never completed |

## Key observations so far

1. **skani search cold-start vs warm.** The first query pays 1308.8 s, of
   which 1300.4 s is `Loading markers` (reading the sketch DB from lustre);
   subsequent queries on the same node drop to 7–45 s. skani search per
   query (warm) is the fastest configuration measured.
2. **FastANI cannot complete one-to-all against the 65,703-genome database**
   on this hardware. Memory grows roughly linearly with thread count
   (~38 GB per thread): killed at 128 GB (t=16), 300 GB (t=8) and 500 GB
   (t=32). Even a single-threaded run exceeded the 64 GB node limit
   (peak RSS 65.4 GB, OOM after ~1.5 h). FastANI is therefore not usable
   for one-to-all searches of this scale under our cluster constraints.
3. **syn2bani's from-FASTA dist is digestion-bound** (~4.1 h per query for
   65,703 references; the 4-enzyme digestion is the cost). Its search mode
   reloads and re-indexes all 65,703 sketches per invocation (~1.5–2.3 h
   cold, ~1.4 h warm) — an engineering gap versus skani's packed marker
   index (persistent index / packed DB is future work). The SV metrics come
   at no extra cost in either mode.
4. **syn2bani search can report tiny shared islands as high-ANI hits.** Two
   of 65,703 hits for the E. coli query showed ani_gated = 99.9999 with
   af_query = 0.0001 (a few hundred perfectly matching bp; skani reports
   nothing for those pairs). Ranking hits by ANI therefore requires an
   aligned-fraction floor (af >= 0.5 used for the top-hit extraction). With
   the floor, syn2bani and skani agree on the top non-self hit for 4/7
   queries. The remaining 3/7 differ at the strain/accession level but are
   taxonomically close: *H. pylori* top hits belong to the same species
   (different clades), while *S. aureus* and *S. rimosus* top hits fall in
   the same genus (`gtdb_one_to_all_tophits.tsv`). A coverage gate in
   `syn2bani search` output is filed as a tool improvement.
5. **SV stage on the top non-self hit** (`rows_sv`): syn2bani struct takes
   0.10–0.54 s and <=15 MB per pair; dnadiff takes 6.1–38.9 s and 79–327 MB
   on the same pairs. 60–90x faster, and struct reuses the digestion that
   the ANI estimate already computed.

## Tool versions / environment

skani 0.3.2, syn2bani 0.1.0 (commit 0b58795), FastANI (conda env fastani),
DNAdiff 1.3. HPC intel partition (16 CPUs, 64 GB) unless noted; fastANI
retries on condo_amd (up to 500 GB).
