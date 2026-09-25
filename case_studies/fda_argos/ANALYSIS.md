# Structural discordance within ANI-identical FDA-ARGOS *E. coli* pairs

Analysis date: 2026-09-25. Pair set: 103 pairs with skani **0.3.2** ANI ≥ 99.9
(`pairs99.9_ecoli.tsv`, regenerated from `skani_triangle_ecoli_v032.tsv`).
Tool: local `syn2bani struct` v0.1.0 (`/Users/macstudio/Downloads/Syn2bANI/target/release/syn2bani`),
per pair: `syn2bani struct <query>.fna <ref>.fna --bed --circular <all ref contigs> -o struct_pairs99.9/<ref>__<query>.bed`
(reference = alphabetically first accession; query = second). All 103 runs exit 0.

## Headline

- **101 of 103 ANI ≥ 99.9 pairs (98%) carry ≥ 1 structural-variant call.**
- Total events: 2,718 — INV 909, INS 707, DEL 661, TRA 441.
- Per pair: median 25 SVs (range 0–61); median summed affected reference span 1.44 Mbp
  (range 0–8.1 Mbp; events can overlap, so sums are an upper bound of unique bp).
- Per-event span: median 4.1 kb, p90 255 kb, max 2.29 Mb; 384 events ≥ 100 kb
  (340 INV + 44 TRA). **83/103 pairs have ≥ 1 event ≥ 100 kb.**
- Even at ANI ≥ 99.99 (24 pairs), 22 carry ≥ 1 SV.
- Not driven by assembly fragmentation alone: the 32 pairs where *both* genomes are
  complete/chromosome-level still show median 15 SVs and 1.44 Mbp affected
  (contig-involved pairs: median 32 SVs / 1.49 Mbp).
- Only 2 pairs are structurally silent — both are ANI 100.0 pairs from the
  GCA_01904xxxxx re-assembly batch (FDAARGOS_1380/1381, FDAARGOS_1378/1379).

## Five most extreme pairs (by summed affected bp)

| pair (ref vs query) | ANI | n SV | affected bp | event mix |
|---|---|---|---|---|
| GCF_016903335.1 vs GCF_016904195.1 | 99.93 | 29 | 8,132,483 | INV-heavy |
| GCF_016903835.1 vs GCF_016904075.1 | 99.90 | 28 | 6,469,239 | INV-heavy |
| GCF_016904195.1 vs GCF_016904535.1 | 99.92 | 31 | 6,233,958 | INV-heavy |
| GCF_002206405.2 vs GCF_017081115.1 | 99.93 | 58 | 5,783,348 | DEL:10 INS:21 INV:23 TRA:4 |
| GCF_002208865.2 vs GCF_017081335.1 | 99.91 | 61 | 5,512,474 | DEL:17 INS:19 INV:21 TRA:4 |

## Metadata concordance (from `metadata_ecoli.tsv`)

- The three most extreme pairs sit in the USA:WA clinical cluster (FDAARGOS_1276–1332)
  and their AMR genotype differences match the structural calls: e.g.
  GCF_016903335.1 vs GCF_016904195.1 differ in `aadA2, dfrA12, mph(A), mrx(A),
  ompC_Q76Ter, sul1` — mobile-element/cassette content consistent with the observed
  INS/DEL events. GCF_016903835.1 vs GCF_016904075.1 has identical AMR genotype
  despite 6.5 Mbp of rearranged span (pure INV differences).
- The FDAARGOS_292/293 vs FDAARGOS_1238/1351/1355 cross-cluster pairs differ in
  `catA1`/`emrD` carriage and geography (USA:DC vs USA:MD).
- Take-home for the paper: within-clone ANI ≥ 99.9 (even ≥ 99.99) does **not** imply
  structural identity — large inversions and AMR-cassette gains/losses are routine,
  which is exactly the gap Syn2bANI's struct track is designed to fill.

## Caveats

- `total_affected_bp` sums per-event reference spans (BED end−start); overlapping
  events double-count shared bases (max-pair sum 8.1 Mb > 5 Mb genome → overlap exists).
- TRA calls and some large INV calls in pairs involving contig-level assemblies may
  include contig-ordering artifacts; the `--circular` origin filter was applied
  (2–3 artifacts filtered per pair per the run logs), but fragmented-assembly
  artifacts are not fully excluded. The both-complete subset (32 pairs) is the
  conservative core: 31/32 have ≥ 1 SV.
- BED score column = anchor support (left+right), not event size; see
  `struct_summary.tsv` and per-pair BEDs in `struct_pairs99.9/`.
