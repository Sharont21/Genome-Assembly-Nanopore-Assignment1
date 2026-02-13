#!/bin/bash

module load StdEnv/2023
module load apptainer/1.4.5

apptainer exec docker://ontresearch/medaka:latest \
  medaka_consensus \
    -i filtered/SRR32410565.filtered.fastq.gz \
    -d assembly/flye/assembly.fasta \
    -o medaka \
    -t 8 \
    -m r1041_e82_400bps_sup_v5.2.0
