#!/bin/bash
# Pileup-based variant calling using bcftools

module load StdEnv/2023
module load samtools
module load bcftools

bcftools mpileup \
  -f reference/GCF_000006945.2_ASM694v2_genomic.fna \
  -Ou alignment/reads/reads_vs_ref.sorted.bam \
| bcftools call \
  -mv \
  -Oz \
  -o bcftools_variants.vcf.gz

bcftools index bcftools_variants.vcf.gz
