#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------
# Download genomes
# ----------------------------------------

# COVID genome
curl -L \
  "https://figshare.com/ndownloader/files/22700666" \
  -o covid.fasta

# E. coli genome
curl -L \
  "https://www.ebi.ac.uk/ena/browser/api/fasta/U00096.3?download=true" \
  -o ecoli.fasta

# ----------------------------------------
# Generate random 1,000,000 bp DNA sequence
# ----------------------------------------

python3 - <<'PY'
import random

length = 1_000_000
bases = "ACGT"

with open("random.fasta", "w") as f:
    f.write(">random_sequence_1Mbp\n")

    seq = ''.join(random.choices(bases, k=length))

    # wrap lines at 80 chars
    for i in range(0, length, 80):
        f.write(seq[i:i+80] + "\n")
PY

echo "Done:"
echo "  covid.fasta"
echo "  ecoli.fasta"
echo "  random.fasta"