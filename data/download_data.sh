#!/bin/bash

curl -L -o ncbi_dataset.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_009914755.1/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"

unzip ncbi_dataset.zip 

cp ncbi_dataset/data/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna ./human_genome.fasta

#!/bin/bash

input="human_genome.fasta"
outdir="chromosomes"

mkdir -p "$outdir"

awk -v outdir="$outdir" '
BEGIN {
    current = ""
}

/^>/ {
    current = ""

    # Match chromosome names 1-22, X, or Y
    if (match($0, /chromosome ([0-9]+|X|Y)/, arr)) {
        chr = tolower(arr[1])
        current = outdir "/ch" chr ".fa"
    }

    # Skip non-chromosome entries
    nextfile_flag = (current == "")
}

{
    if (!nextfile_flag) {
        print >> current
    }
}
' "$input"