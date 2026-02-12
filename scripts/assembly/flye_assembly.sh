#!/bin/bash

apptainer exec ~/containers/flye.sif \
  flye \
    --nano-hq filtered/SRR32410565.filtered.fastq.gz \
    --genome-size 4.8m \
    --out-dir assembly/flye \
    --threads 8

