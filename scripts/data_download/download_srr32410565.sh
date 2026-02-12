#!/bin/bash

# Download raw ONT data from NCBI SRA
# Accession: SRR32410565
# Run on Nibi (Digital Research Alliance of Canada)

module load StdEnv/2023
module load sra-toolkit

mkdir -p data
cd data

fasterq-dump SRR32410565 --threads 8
gzip SRR32410565.fastq
