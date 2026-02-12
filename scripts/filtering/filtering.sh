apptainer exec ~/containers/seqkit.sif \
  seqkit seq \
    --min-len 1000 \
    --min-qual 10 \
    data/SRR32410565.fastq.gz \
  > filtered/SRR32410565.filtered.fastq
