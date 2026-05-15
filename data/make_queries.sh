#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------
# For every .fasta file:
#   generate:
#     <name>20queries.txt
#     <name>50queries.txt
#
# Also generate:
#   randomqueries.txt
# ----------------------------------------

python3 - <<'PY'
import random
import glob
import os

NUM_QUERIES = 1000

def read_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line.upper())
    return ''.join(seq)

def sample_queries(seq, length, count):
    if len(seq) < length:
        raise ValueError(f"Sequence shorter than query length {length}")

    out = []
    max_start = len(seq) - length

    for _ in range(count):
        start = random.randint(0, max_start)
        out.append(seq[start:start+length])

    return out

# ----------------------------------------
# Process FASTA files
# ----------------------------------------

for fasta in glob.glob("*.fasta"):
    basename = os.path.splitext(fasta)[0]

    print(f"Processing {fasta}...")

    seq = read_fasta(fasta)

    q20 = sample_queries(seq, 20, NUM_QUERIES)
    q50 = sample_queries(seq, 50, NUM_QUERIES)

    with open(f"{basename}20queries.txt", "w") as f:
        f.write("\n".join(q20) + "\n")

    with open(f"{basename}50queries.txt", "w") as f:
        f.write("\n".join(q50) + "\n")

# ----------------------------------------
# Generate fully random queries
# ----------------------------------------

bases = "ACGT"

with open("randomqueries.txt", "w") as f:
    for _ in range(NUM_QUERIES):
        q = ''.join(random.choices(bases, k=10))
        f.write(q + "\n")

print("Done.")
PY