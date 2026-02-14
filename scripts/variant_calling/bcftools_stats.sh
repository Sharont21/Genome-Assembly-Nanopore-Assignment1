#!/bin/bash
module load StdEnv/2023
module load bcftools

bcftools view -v snps clair3/merge_output.vcf.gz | wc -l
bcftools view -v indels clair3/merge_output.vcf.gz | wc -l
