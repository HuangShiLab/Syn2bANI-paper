# HPC layout for Syn2bANI-paper

The canonical HPC checkout is:

```text
/lustre1/g/aos_shihuang/Syn2bANI-paper
```

This checkout tracks `origin/main`. Do not create additional `Syn2bANI-*`
working copies at the group-directory root.

HPC-only inputs, benchmark scratch space and census working data are stored in
the ignored directory:

```text
/lustre1/g/aos_shihuang/Syn2bANI-paper/hpc_workspace/
```

Current subdirectories are:

| HPC path under `hpc_workspace/` | Contents |
|---|---|
| `code/Syn2bANI` | HPC-local Syn2bANI code/build workspace |
| `benchmarks/Syn2bANI-bench` | Legacy Syn2bANI benchmark clone and artifacts |
| `benchmarks/Syn2bANI-hi95` | Legacy 95–99.5% ANIm expansion workspace |
| `benchmarks/Syn2bANI-paper-bench` | GTDB/SynTracker benchmark scratch space |
| `benchmarks/syn2b_benchmarks` | Legacy Syn2b benchmark workspace |
| `benchmarks/validation_mid_ani` | Mid-ANI validation workspace |
| `census/syn2b_census` | GTDB within-species census working data |
| `census/syntracker_validation` | SynTracker structural-validation data |

Build caches, SLURM logs and other generated scratch files are removed from
these workspaces as needed. Census outputs and final summary tables are retained
until their Zenodo archive is complete.
