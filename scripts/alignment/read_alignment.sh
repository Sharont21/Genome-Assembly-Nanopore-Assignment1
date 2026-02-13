#!/bin/bash

# Read-to-reference alignment for variant calling

cd /scratch/sharont/salmonella_nanopore || exit 1

module load StdEnv/2023
module load minimap2
module load samtools

mkdir -p alignment/reads

minimap2 -a -x map-ont \
  reference/GCF_000006945.2_ASM694v2_genomic.fna \
  filtered/SRR32410565.filtered.fastq.gz \
  | samtools view -bS \
  | samtools sort -o alignment/reads/reads_vs_ref.sorted.bam

samtools index alignment/reads/reads_vs_ref.sorted.bam
