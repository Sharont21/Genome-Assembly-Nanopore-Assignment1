# Genome-Assembly-Nanopore-Assignment1
Genome assembly and reference comparison of _Salmonella enterica_ using Oxford Nanopore long-read sequencing

# Introduction
_De novo_ genome assembly is used to reconstruct an organism’s genome by computationally stitching together millions of fragmented sequencing reads based on sequence overlap, a process complicated by sequencing errors, uneven coverage, and missing data (Baker, 2012). To overcome these challenges, assemblers rely on multiple overlapping reads covering each genomic region to produce an accurate consensus sequence (Baker, 2012). Recent developments in high-throughput and long-read sequencing platforms have greatly broadened the ability to characterize microbial genomes, allowing researchers to examine genomic structure, diversity, and evolutionary patterns with increased resolution across biological and medical studies (Zhou and Faust, 2025).

Assembling a genome from raw sequencing reads is computationally and biologically challenging. Sequencing reads contain errors that must be distinguished from true biological variation, and assembly algorithms must accurately identify overlaps between reads to reconstruct the original genome sequence (Altermann et al., 2022). In addition, uneven sequencing coverage can result in gaps or reduced contiguity, requiring tradeoffs between assembly completeness and accuracy (Altermann et al., 2022). These challenges are strongly influenced by the sequencing technology used and the error profiles associated with different platforms.

Long-read sequencing technologies such as Oxford Nanopore Technologies (ONT) improve genome assembly contiguity by producing reads that span repetitive regions and structural variants (Bogaerts et al., 2025). Recent advances in ONT R10 (“Q20”) chemistry have reduced base-calling error rates, enabling near-complete bacterial genome assemblies in some cases without additional short-read data (Bogaerts et al., 2025). ONT sequencing infers nucleotide sequences from changes in electrical current as DNA passes through a nanopore, with machine-learning-based basecalling translating these signals into long reads that resolve complex genomic regions (Bogaerts et al., 2025). Despite these improvements, ONT data remain more error-prone than short-read sequencing, with residual insertion–deletion errors, particularly in homopolymer regions (Wick and Holt, 2022). Although hybrid polishing approaches using short-read data can improve consensus accuracy, alignment ambiguity in repetitive regions often leaves some errors unresolved, underscoring the need to balance assembly contiguity, accuracy, and computational complexity when selecting an assembly strategy (Wick and Holt, 2022).

Following assembly, alignment to a reference genome provides a framework for assessing assembly quality and identifying genomic variation. Reference-based comparison enables efficient detection of single nucleotide polymorphisms and small insertions or deletions when a closely related reference genome is available (Vo et al., 2013). However, this approach may introduce reference bias, as divergent or novel genomic regions may align poorly or be excluded, potentially obscuring structural variation or strain-specific content.

In this analysis, raw ONT R10 sequencing reads from Salmonella enterica will be assembled into a consensus genome and compared to a reference genome obtained from NCBI. Several assembly and analysis tools were evaluated based on their suitability for long-read bacterial genome assembly, computational efficiency, and accuracy. The selected workflow balances contiguity, error correction, and interpretability, and includes genome assembly, quality assessment, reference alignment, variant calling, and visualization to evaluate both assembly performance and genomic variation.

Several long-read assemblers were evaluated for assembling the Salmonella enterica genome from Oxford Nanopore data, including Flye, Canu, and NECAT. Flye was selected as the primary assembler due to its robustness and consistent performance for long-read de novo bacterial genome assembly (Hussain et al., 2023) and (Wick and Holt, 2019). Flye constructs an assembly graph by connecting error-prone disjoint genomic segments (“disjointigs”) into a unified structure, enabling accurate reconstruction of complex genomic regions, and has been shown to perform particularly well for plasmid assembly (Hussain et al., 2023). While Canu applies extensive read correction and trimming to reduce base-calling errors and can generate reliable assemblies, it has substantially longer runtimes and higher computational demands, making it less efficient for this analysis (Wick and Holt, 2019) and (Hussain et al., 2023). NECAT employs a two-stage correction and assembly strategy designed to address Nanopore-specific errors and can generate high-quality assemblies more rapidly than Canu; however, it often requires additional polishing steps to achieve comparable accuracy (Hussain et al., 2023). Based on robustness, accuracy, and computational efficiency, Flye was chosen as the most suitable assembler for this workflow.

For reference alignment, minimap2 was chosen due to its widespread use and optimized performance for mapping long-read assemblies to closely related reference genomes (Li, 2021). Variant detection and file processing will be conducted using SAMtools and BCFtools for manipulating alignment files and identifying genomic variation. Finally, the Integrative Genomics Viewer (IGV) will be used for visualization, allowing manual inspection of alignments and variants to assess assembly quality and validate detected sequence differences. Together, these tools form a workflow that balances computational efficiency, accuracy, and interpretability for long-read bacterial genome assembly and comparison.

# Methods
##Genome assembly and polishing

Oxford Nanopore R10 sequencing reads (FASTQ format) will be assembled _de novo_ using Flye (v2.9.6), a long-read assembler designed for error-prone Nanopore data. Assembly was performed in high-accuracy Nanopore mode `--nano-hq`, with an estimated genome size of 4.8 Mb `--genome-size 4.8m`, consistent with _Salmonella enterica_ (NBCI). Default parameters will used unless otherwise specified. Flye’s internal polishing steps will applied to improve consensus accuracy prior to downstream analysis.

`flye --nano-hq reads.fastq --genome-size 4.8m --out-dir flye_output --threads 8`


##Reference genome alignment
To evaluate the assembled genome and identify sequence differences, the draft assembly was aligned to a Salmonella enterica reference genome downloaded from the NCBI RefSeq database. Alignment was performed using minimap2 (v2.26), which is well suited for long-read and assembly-to-reference alignments due to its speed and accuracy when handling large, highly similar sequences.

##Alignment processing and variant calling
Alignment files were converted, sorted, and indexed using SAMtools (v1.19.2) to prepare them for variant analysis. Variants relative to the reference genome, including single nucleotide polymorphisms (SNPs) and small insertions and deletions (indels), were identified using BCFtools (v1.19). This reference-based variant calling approach enabled systematic detection of residual errors and true biological variation remaining in the assembled genome.

##Visualization
To visually inspect alignment quality and genomic differences between the assembled genome and the reference, alignments and variant calls were examined using Integrative Genomics Viewer (IGV, v2.19.7). IGV allowed for manual assessment of variant support, coverage consistency, and error-prone regions such as homopolymers.

# References
[1] Zhou, X., & Faust, K. (2025). A high-throughput and time-efficient Nanopore full-length 16S rRNA gene sequencing protocol for synthetic microbial communities. Methods (San Diego, Calif.), 240, 14–20. https://doi.org/10.1016/j.ymeth.2025.04.003

[2] Bogaerts, B., Maex, M., Commans, F., Goeders, N., Van den Bossche, A., De Keersmaecker, S. C. J., Roosens, N. H. C., Ceyssens, P.-J., Mattheus, W., & Vanneste, K. (2025). Oxford Nanopore Technologies R10 sequencing enables accurate cgMLST-based bacterial outbreak investigation of Neisseria meningitidis and Salmonella enterica when accounting for methylation-related errors. Journal of Clinical Microbiology, 63(10), e0041025. https://doi.org/10.1128/jcm.00410-25

[3] Altermann, E., Tegetmeyer, H. E., & Chanyi, R. M. (2022). The evolution of bacterial genome assemblies—Where do we need to go next. Microbiome Research Reports, 1(2), 15. https://doi.org/10.20517/mrr.2022.02

[4] Wick, R. R., & Holt, K. E. (2022). Polypolish: Short-read polishing of long-read bacterial genome assemblies. PLoS Computational Biology, 18(1), e1009802. https://doi.org/10.1371/journal.pcbi.1009802

[5] Zhou, X., & Faust, K. (2025). A high-throughput and time-efficient Nanopore full-length 16S rRNA gene sequencing protocol for synthetic microbial communities. Methods (San Diego, Calif.), 240, 14–20. https://doi.org/10.1016/j.ymeth.2025.04.003

[6] Vo, N. S., Tran, Q., Niraula, N., & Phan, V. (2014). RandAL: a randomized approach to aligning DNA sequences to reference genomes. BMC Genomics, 15(Suppl 5), Article S2. https://doi.org/10.1186/1471-2164-15-S5-S2

[7] Wick, R. R., & Holt, K. E. (2019). Benchmarking of long-read assemblers for prokaryote whole genome sequencing [version 1; peer review: 4 approved]. F1000 Research, 8, 2138. https://doi.org/10.12688/f1000research.21782.1

[8] Safar, H. A., Alatar, F., Nasser, K., Al-Ajmi, R., Alfouzan, W., & Mustafa, A. S. (2023). The impact of applying various de novo assembly and correction tools on the identification of genome characterization, drug resistance, and virulence factors of clinical isolates using ONT sequencing. BMC Biotechnology, 23(1), Article 26. https://doi.org/10.1186/s12896-023-00797-3

[9] Baker, M. (2012). De novo genome assembly: what every biologist should know. Nature Methods, 9(4), 333–337. https://doi.org/10.1038/nmeth.1935

[10] Li, H. (2021). New strategies to improve minimap2 alignment accuracy. Bioinformatics, 37(23), 4572–4574. https://doi.org/10.1093/bioinformatics/btab705

