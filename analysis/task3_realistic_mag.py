#!/usr/bin/env python3
"""Task 3: Realistic MAG validation with contamination, chimerism, duplication, misassembly."""
import random, subprocess, csv, os
from pathlib import Path

random.seed(42)

BENCH = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI_benchmark_realistic")
BENCH.mkdir(exist_ok=True)
ECOLI = Path("/Users/shihuang/test_data/GCF_002953055.1_ASM295305v1_genomic.fna")
SYN2B = Path("/Users/shihuang/Documents/kimi/workspace/Syn2bANI/target/release/syn2bani")
REF = BENCH / "reference.fasta"

def parse_fasta(path):
    seqs = {}
    cid = None; buf = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if cid: seqs[cid] = ''.join(buf)
                cid = line[1:].split()[0]; buf = []
            elif line: buf.append(line)
        if cid: seqs[cid] = ''.join(buf)
    return seqs

def write_fasta(path, seqs):
    with open(path, 'w') as f:
        for sid, seq in seqs.items():
            f.write(f'>{sid}\n')
            for i in range(0, len(seq), 80): f.write(seq[i:i+80]+'\n')

ref_seqs = parse_fasta(ECOLI)
ref_seq = list(ref_seqs.values())[0]
write_fasta(REF, {"reference": ref_seq})

print("="*60)
print("Task 3: Realistic MAG Scenarios")
print("="*60)

scenarios = {}

# 1. Contamination (add random DNA as foreign species)
for contam_frac in [0.0, 0.05, 0.10, 0.20]:
    contam_len = int(len(ref_seq) * contam_frac)
    contam = ''.join(random.choices('ACGT', k=contam_len))
    # Insert contamination blocks at random positions
    modified = ref_seq
    if contam_len > 0:
        pos = random.randint(0, len(modified) - contam_len)
        modified = modified[:pos] + contam + modified[pos:]
    name = f"contam_{contam_frac:.2f}"
    write_fasta(BENCH / f"{name}.fasta", {name: modified})
    scenarios[name] = {"type": "contamination", "param": contam_frac, "desc": f"{contam_frac*100:.0f}% random DNA contamination"}
    print(f"  {name}: {contam_len} bp contamination ({contam_frac*100:.0f}%)")

# 2. Chimerism (misorder contigs from different positions)
def create_chimeric(seq, n_breaks=5):
    """Break genome into pieces and shuffle order."""
    breaks = sorted(random.sample(range(1000, len(seq)-1000), n_breaks))
    pieces = []
    prev = 0
    for b in breaks:
        pieces.append(seq[prev:b])
        prev = b
    pieces.append(seq[prev:])
    random.shuffle(pieces)
    return ''.join(pieces), n_breaks

for n_breaks in [1, 3, 5, 10, 20]:
    chim_seq, _ = create_chimeric(ref_seq, n_breaks)
    name = f"chim_{n_breaks}"
    write_fasta(BENCH / f"{name}.fasta", {name: chim_seq})
    scenarios[name] = {"type": "chimerism", "param": n_breaks, "desc": f"{n_breaks} breakpoints, shuffled"}
    print(f"  {name}: {n_breaks} breakpoints, shuffled")

# 3. Duplication (copy some segments)
def create_duplicated(seq, dup_frac=0.05):
    """Duplicate random segments."""
    modified = seq
    dup_len = int(len(seq) * dup_frac)
    while dup_len > 1000:
        seg_len = min(dup_len, random.randint(1000, 10000))
        start = random.randint(0, len(modified) - seg_len)
        segment = modified[start:start+seg_len]
        insert_pos = random.randint(0, len(modified))
        modified = modified[:insert_pos] + segment + modified[insert_pos:]
        dup_len -= seg_len
    return modified

for dup_frac in [0.0, 0.05, 0.10, 0.20]:
    dup_seq = create_duplicated(ref_seq, dup_frac)
    name = f"dup_{dup_frac:.2f}"
    write_fasta(BENCH / f"{name}.fasta", {name: dup_seq})
    scenarios[name] = {"type": "duplication", "param": dup_frac, "desc": f"{dup_frac*100:.0f}% sequence duplication"}
    print(f"  {name}: {dup_frac*100:.0f}% duplication")

# 4. Assembly error (introduce SNPs/indels into contigs)
def introduce_errors(seq, error_rate=0.001):
    """Introduce SNPs and small indels."""
    seq_list = list(seq)
    for i in range(len(seq_list)):
        if random.random() < error_rate:
            old = seq_list[i]
            seq_list[i] = random.choice([b for b in 'ACGT' if b != old])
    # Small indels
    result = []
    i = 0
    while i < len(seq_list):
        if random.random() < error_rate * 0.1:
            if random.random() < 0.5:
                i += 1; continue
            else:
                result.append(random.choice('ACGT'))
        result.append(seq_list[i])
        i += 1
    return ''.join(result)

for err_rate in [0.0, 0.0001, 0.0005, 0.001, 0.002]:
    err_seq = introduce_errors(ref_seq, err_rate)
    name = f"err_{err_rate:.4f}"
    write_fasta(BENCH / f"{name}.fasta", {name: err_seq})
    scenarios[name] = {"type": "assembly_error", "param": err_rate, "desc": f"{err_rate*100:.2f}% SNP/indel error"}
    print(f"  {name}: {err_rate*100:.2f}% assembly error")

# 5. Combined realistic MAG (all defects at once)
realistic_seq = ref_seq
realistic_seq, _ = create_chimeric(realistic_seq, 10)
realistic_seq = introduce_errors(realistic_seq, 0.001)
realistic_seq = create_duplicated(realistic_seq, 0.05)
# Add contamination
contam_len = int(len(realistic_seq) * 0.05)
contam = ''.join(random.choices('ACGT', k=contam_len))
pos = random.randint(0, len(realistic_seq) - contam_len)
realistic_seq = realistic_seq[:pos] + contam + realistic_seq[pos:]
# Fragment into realistic MAG contigs
mean_frag = 10000
frags = []
pos = 0
while pos < len(realistic_seq):
    fl = max(500, int(random.expovariate(1.0/mean_frag)))
    if pos + fl > len(realistic_seq): fl = len(realistic_seq) - pos
    frags.append(realistic_seq[pos:pos+fl])
    pos += fl

frag_seqs = {f"contig_{i}": f for i, f in enumerate(frags)}
write_fasta(BENCH / "realistic_mag.fasta", frag_seqs)
scenarios["realistic_mag"] = {"type": "combined", "param": 0, "desc": "Realistic MAG with all defects"}
print(f"  realistic_mag: {len(frags)} contigs, all defects combined")

# Run Syn2bANI on all scenarios
print("\n" + "="*60)
print("Running Syn2bANI on realistic MAG scenarios")
print("="*60)

results = []
for name in sorted(scenarios.keys()):
    q_path = BENCH / f"{name}.fasta"
    cmd = [str(SYN2B), "dist", str(q_path), str(REF)]
    try:
        output = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        for line in output.stdout.strip().split('\n'):
            if line.startswith('/'):
                parts = line.split('\t')
                if len(parts) >= 9:
                    ani = float(parts[2])
                    af_q = float(parts[3])
                    af_r = float(parts[4])
                    shared = int(parts[7])
                    results.append({
                        'name': name,
                        'scenario': scenarios[name]['type'],
                        'param': scenarios[name]['param'],
                        'ani': ani,
                        'af_q': af_q,
                        'af_r': af_r,
                        'shared_tags': shared,
                        'desc': scenarios[name]['desc']
                    })
                    print(f"  {name}: ANI={ani:.4f}, AF={af_q:.4f}, tags={shared}")
    except Exception as e:
        print(f"  ERROR {name}: {e}")

# Save results
with open(BENCH / "realistic_mag_results.csv", 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['name', 'scenario', 'param', 'ani', 'af_q', 'af_r', 'shared_tags'])
    for r in results:
        writer.writerow([r['name'], r['scenario'], r['param'], r['ani'], r['af_q'], r['af_r'], r['shared_tags']])

print(f"\nResults saved to {BENCH / 'realistic_mag_results.csv'}")
print(f"Total scenarios tested: {len(results)}")
