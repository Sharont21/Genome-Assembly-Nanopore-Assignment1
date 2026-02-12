apptainer exec ~/containers/nanoplot.sif \
  NanoPlot \
    --fastq data/SRR32410565.fastq.gz \
    --outdir qc_before \
    --threads 4 \
    --loglength
