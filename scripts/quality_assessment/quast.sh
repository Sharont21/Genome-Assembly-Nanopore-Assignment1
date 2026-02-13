#!/bin/bash

conda activate quast_env

quast.py \
  assembly/flye/assembly.fasta \
  -r reference/GCF_000006945.2_ASM694v2_genomic.fna \
  -o quast \
  --threads 4

conda deactivate
