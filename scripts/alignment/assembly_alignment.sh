#!/bin/bash

module load StdEnv/2023
module load minimap2
module load samtools

minimap2 -ax asm5 \
  reference/GCF_000006945.2_ASM694v2_genomic.fna \
  medaka/consensus.fasta \
  > alignment/assembly/assembly_vs_ref.sam

samtools view -bS alignment/assembly/assembly_vs_ref.sam \
  | samtools sort -o alignment/assembly/assembly_vs_ref.sorted.bam

samtools index alignment/assembly/assembly_vs_ref.sorted.bam
