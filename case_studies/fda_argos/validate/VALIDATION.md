# Validation of syn2bani struct calls against minimap2 (FDA-ARGOS E. coli, 103 pairs)

Date: 2026-09-25. Method: `minimap2 -cx asm20` (and `--cs` re-run) per pair,
reference = alphabetically first accession (`validate/pafs/`, `validate/pafs_cs/`).
Per-pair classification uses only the longest (chromosome) reference contig:
minus_fraction = aligned ref bases on − strand / (+ + −); class =
`global_flip_or_rotation` (>0.8), `partial_inversion` (0.2–0.8), `collinear_local` (≤0.2).
Machine-readable per-pair results: `validation_summary.tsv` (script: `validate/analyze_pafs.py`).

## Why: the flip/rotation artifact

Manual check of the most extreme pair (FDAARGOS_401 vs 772, ANI 100) showed the query
chromosome aligns entirely on the minus strand in exactly two blocks
(0–3,437,929 | 3,437,929–5,311,821) — a reverse-complement + circular-rotation
assembly convention, not a biological 1.87 Mb inversion. The struct INV call sat exactly
at the rotation breakpoint. Our `--circular` filter only removes calls spanning >50% of a
contig, so rotation artifacts presenting below that threshold survived.

## Classification of the 103 pairs

| class | pairs | interpretation |
|---|---|---|
| global_flip_or_rotation | 41 (40%) | assembly-wide strand/origin convention; large INV/TRA calls here are artifacts |
| partial_inversion | 16 (16%) | genuine large inversion(s) (minus-strand blocks on plus background) |
| collinear_local | 46 (45%) | collinear; struct calls are local events |

Odd-structure flags: 14 pairs where the chromosome carries <50% of aligned bases,
2 with no dominant chromosome (draft), 17 with >30 alignment blocks (fragmented).

## Event-level accounting (2,718 struct calls)

| category | events | verdict |
|---|---|---|
| INV/TRA inside flipped background or at rotation junctions | 355 | **artifact** (coordinate convention) |
| indels at 2-block rotation junctions | 3 | artifact |
| chromosomal INS/DEL ≥1 kb verified in PAF cs strings | 876 | **genuine** |
| rearrangements strand-confirmed (minus block, or plus-island in flipped pair) or uncontradicted | 445 | **genuine** (incl. all 16 partial-inversion pairs) |
| off-chromosome (plasmid/contig) content differences | 880 | **genuine** |
| small indels below cs-verification threshold | 159 | likely genuine, unverified |

Deep check of the 10 most extreme flipped pairs: their big INV/TRA calls are all
artifact-zone (rotation junctions), but every one still carries 4–32 cs-verified indels
and/or plasmid-content differences (e.g. GCF_002208865.2 vs GCF_017081375.1: 19 artifact
rearrangement calls vs 32 cs-verified indels).

## Corrected headline numbers

- **101/103 ANI ≥ 99.9 pairs (98%) carry ≥ 1 genuine structural difference** — the
  pair-level headline from `ANALYSIS.md` survives validation.
- **What changes:** 358/2,718 struct calls (13%) are orientation/origin artifacts, and
  these were the *largest* calls — the "1.87–8.1 Mb affected" claims for flipped pairs
  (incl. the top-5 extreme pairs in ANALYSIS.md) are convention artifacts. Genuine
  per-pair structural burden is dominated by indels (median ~4 kb, up to 52 kb
  prophage/cassette-scale), plasmid-content differences, and — in the 16
  partial-inversion pairs — real large inversions.
- 41/41 flipped pairs still differ by genuine indels/plasmid content; flip status itself
  is annotation, not biology.
- 2/103 pairs are fully structurally silent (both ANI-100.0 re-assembly pairs from the
  GCA_01904xxxxx batch).

## Flagship verified pairs (ANI ≥ 99.9, both complete genomes, PAF-verified)

| pair (strains) | ANI | verified content | interpretation |
|---|---|---|---|
| GCF_016889265.1 vs GCF_016890045.1 (FDAARGOS_1264/1265) | 99.93 | 22 cs-verified indels (51/45/40/38 kb) + plasmid events; AMR diff: **blaCTX-M-15, blaOXA-1, aac(6')-Ib-cr5, catB3, tet(A)** | ESBL cassette + prophage differences within a clone |
| GCF_002209105.2 vs GCF_002393365.1 (FDAARGOS_348/403) | 99.97 | 16 cs-verified indels (16 kb, 15 kb); AMR diff: **blaCTX-M-15** | CTX-M-15 acquisition in O104:H4-lineage pair |
| GCF_002209105.2 vs GCF_006364695.1 (FDAARGOS_348/772) | 100.0 | 48 kb INS (prophage/cassette), 16 kb INS, 9 kb DEL; AMR diff: **tet(A), tet(K)** | ANI-identical pair differing by tet cassette + prophage |
| GCF_016903335.1 vs GCF_016904535.1 (FDAARGOS_1295/1297) | 99.93 | 46 kb DEL + 46 kb INS (prophage swap); AMR diff: **aadA5, aph(3'')-Ib, aph(6)-Id, blaTEM-1, dfrA17, sul1, sul2** | integron/cassette exchange, USA:WA clinical cluster |
| GCF_016889045.1 vs GCF_016904055.1 (FDAARGOS_1255/1293) | 99.95 | 4 prophage-scale indels (52/50/46/41 kb); acrR frameshift diff | prophage turnover within clone |
| GCF_002393365.1 vs GCF_006364695.1 (FDAARGOS_403/772) | 99.90 | 48 kb INS + 15/10 kb DELs; AMR diff: **blaCTX-M-15, tet(A), tet(K)** | ESBL + tet cassette differences |

(Plus the 16 partial-inversion pairs as a class: genuine minus-strand-confirmed large
inversions at ANI ≥ 99.9.)

## Caveats

- cs-verification threshold: indels ≥1 kb on the chromosome; smaller struct indels
  (159 events) are counted "likely genuine, unverified" rather than confirmed.
- TRA calls involving contig-level assemblies remain partly unresolvable (17 pairs with
  >30 blocks); the both-complete subset is the clean core.
- minus_fraction classes use the reference chromosome only; plasmid strands ignored.
