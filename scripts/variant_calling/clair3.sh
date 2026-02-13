#!/bin/bash

module load StdEnv/2023
module load minimap2
module load samtools
module load apptainer/1.4.5

# Index the reference genome (required for Clair3)
samtools faidx reference/GCF_000006945.2_ASM694v2_genomic.fna

# Run Clair3 variant calling
apptainer exec docker://hkubal/clair3:latest \
  run_clair3.sh \
    --bam_fn alignment/reads/reads_vs_ref.sorted.bam \
    --ref_fn reference/GCF_000006945.2_ASM694v2_genomic.fna \
    --threads 8 \
    --platform ont \
    --model_path /opt/models/r1041_e82_400bps_sup_v500 \
    --output clair3
