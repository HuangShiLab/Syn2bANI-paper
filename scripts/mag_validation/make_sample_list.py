#!/usr/bin/env python3
"""Generate sample_list.tsv: 25 random CAMI2 strain-madness samples (seed 42)
+ all 10 CAMI2 marine samples. Columns: dataset, sample.
(reads_path is derived on the HPC side via common.sh:reads_path.)
"""
import random

random.seed(42)
strain = sorted(random.sample(range(100), 25))
marine = list(range(10))

with open("lists/sample_list.tsv", "w") as f:
    f.write("dataset\tsample\n")
    for n in strain:
        f.write(f"strain\tsample_{n}\n")
    for n in marine:
        f.write(f"marine\tsample_{n}\n")
print("strain samples:", strain)
print("marine samples:", marine)
