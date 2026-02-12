#!/bin/bash

# NanoPlot QC before filtering
# Tool: NanoPlot v1.46.2 (Apptainer container)

apptainer exec ~/containers/nanoplot.sif \
  NanoPlot \
    --fastq data/SRR32410565.fastq.gz \
    --outdir qc_before \
    --threads 4 \
    --loglength
