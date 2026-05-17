# FM-index for RePair compressed genomes

This repository contains a first implementation of the FM-index for RePair-compressed genomes.

## Requirements

Running this requires:

 - gcc 11+
 - Make
 - Linux OS
 - Python 3.11+

## Reproducing results

Clone the github:

```bash
git clone https://github.com/ambrozak/RPFMIndex.git

cd RPFMIndex
```

Compile the binaries in [./src](./src) using make:

```bash
cd src
make
cd ..
```

Then download and prepare the data in in [./data](./data):

```bash
cd data
bash get_data.sh
bash make_queries.sh
```

Lastly, navigate to the desired dataset folder in [./eval](./eval) and automatically run the queries.

```bash
cd ./eval/covid
# Compression amount must be at least 3, and is the bits allowed per symbol (this grows exponentially!)
# bash build_indexes.sh [COMPRESSION_AMOUNT]
bash build_indexes.sh 5
# Make sure you use the same compression amount for both commands! Index files are named based on this parameter.
# bash query_indexes [COMPRESSION_AMOUNT] [TRIALS] 
bash query_indexes 5 3
```
