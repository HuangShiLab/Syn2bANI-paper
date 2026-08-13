# Machine & Benchmark Configuration

## Machine
- Model: Mac Studio (Mac16,9, Z1CD00236ZP/A)
- Chip: Apple M4 Max, 16 cores (12 performance + 4 efficiency)
- RAM: 128 GB (137438953472 bytes, `sysctl hw.memsize`)
- Logical CPUs (`sysctl hw.ncpu`): 16
- OS: macOS 26.5.1 (build 25F80)
- Benchmark date: 2026-08-13 (UTC)

## Tool versions
- syn2bani 0.1.0 — `/Users/macstudio/Downloads/Syn2bANI/target/release/syn2bani`
  - built from Syn2bANI main @ `69ce9f4`
  - default enzyme panel BcgI,AlfI,AloI,FalI (no `-e` passed to `ani`)
- skani 0.1.0 — `/Users/macstudio/.cargo/bin/skani`
- FastANI version 1.33 — `/opt/homebrew/bin/fastani`

## Threads
- syn2bani: `-t 0` (auto → all cores)
- skani: `-t 16`
- FastANI: `-t 16`

## Genomes
- 14 complete Enterobacteriaceae chromosomes from
  `Syn2bANI/prototype/realbench/genomes/*.fasta` (note: directory contains 14
  assemblies, not 13 — all were used)
- 8 real draft E. coli assemblies from
  `Syn2bANI/prototype/draftbench/drafts/GCA_*.fasta` (forward orientation only;
  the `.rc.fasta` reverse-complement duplicates were excluded)
- Largest subset: n = 22 (mixed complete + draft)
- Nested subsets: n = 2, 5, 10, 15, 22 in a fixed order (E. coli completes
  first, then other genera, then drafts); membership in `genome_subsets.tsv`

## Design
- All-vs-all (n x n, including self pairs) per subset; n_pairs = n^2.
- 3 repetitions per (tool, mode, n); all reps kept in the TSVs, medians used
  for summaries.
- Timing: `/usr/bin/time -l` (macOS); wall time = `real` seconds, peak memory =
  `maximum resident set size` (bytes, converted to MiB).
- Tool stdout/stderr redirected to per-run logs under `logs/`; per-run time
  output in `logs/*.log.time`.

## Exact command lines (per subset, list files under `lists/`)
- syn2bani, FASTA mode (re-digests in memory each run):
  `syn2bani ani --ql lists/fasta_n<N>.txt --rl lists/fasta_n<N>.txt -t 0 -o <out.tsv>`
- syn2bani, sketch step (timed separately, 3 reps):
  `syn2bani sketch --enzymes BcgI,AlfI,AloI,FalI -t 0 -o sketches/syn2bani_n<N> <fastas...>`
  (the panel is passed explicitly to `sketch` because its own default is BcgI
  only; the list matches the `ani` default panel exactly)
- syn2bani, sketch-reuse mode:
  `syn2bani ani --ql lists/s2ba_n<N>.txt --rl lists/s2ba_n<N>.txt -t 0 -o <out.tsv>`
- skani, sketch step:
  `skani sketch -t 16 -o sketches/skani_n<N> <fastas...>`
- skani, dist on sketches:
  `skani dist -t 16 --ql lists/skani_n<N>.txt --rl lists/skani_n<N>.txt -o <out.tsv>`
- FastANI (no sketch concept, run from FASTA):
  `fastani --ql lists/fasta_n<N>.txt --rl lists/fasta_n<N>.txt -t 16 -o <out.tsv>`
  (all FastANI defaults: k=16, fragLen=3000; not capped — n=22 completed well
  under the 15 min limit)

## Sketch artifact sizes
- `du -sk` totals per subset in `sketch_benchmark.tsv` (`total_size_kb`),
  recorded for each rep; artifacts under `sketches/syn2bani_n<N>/` (one `.s2ba`
  per genome) and `sketches/skani_n<N>/` (one `.sketch` per genome +
  `markers.bin`).

## Reproduce
`bash run_benchmark.sh` (full driver with the fixed genome order embedded).
