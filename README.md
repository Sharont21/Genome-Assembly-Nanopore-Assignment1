# Genome-Assembly-Nanopore-Assignment1
Genome assembly and reference comparison of _Salmonella enterica_ using Oxford Nanopore long-read sequencing

# Introduction
_De novo_ genome assembly is used to reconstruct an organism’s genome by computationally stitching together millions of fragmented sequencing reads based on sequence overlap, a process complicated by sequencing errors, uneven coverage, and missing data [1]. To overcome these challenges, assemblers rely on multiple overlapping reads covering each genomic region to produce an accurate consensus sequence [1]. Recent developments in high-throughput and long-read sequencing platforms have greatly broadened the ability to characterize microbial genomes, allowing researchers to examine genomic structure and diversity with increased resolution across biological and medical studies [2].

Assembling a genome from raw sequencing reads is computationally and biologically challenging. Sequencing reads contain errors that must be distinguished from true biological variation, and assembly algorithms must accurately identify overlaps between reads to reconstruct the original genome sequence [3]. In addition, uneven sequencing coverage can result in gaps or reduced contiguity, requiring tradeoffs between assembly completeness and accuracy [3]. These challenges are strongly influenced by the sequencing technology used and the error profiles associated with different platforms.

Long-read sequencing technologies such as Oxford Nanopore Technologies (ONT) improve genome assembly contiguity by producing reads that span repetitive regions and structural variants [4]. Recent advances in ONT R10 (“Q20”) chemistry have reduced base-calling error rates, enabling near-complete bacterial genome assemblies in some cases without additional short-read data [4]. ONT sequencing infers nucleotide sequences from changes in electrical current as DNA passes through a nanopore, with machine-learning-based basecalling translating these signals into long reads that resolve complex genomic regions [3,4]. Despite these improvements, ONT data remain more error-prone than short-read sequencing, with residual insertion–deletion errors, particularly in homopolymer regions [4,5]. Although hybrid polishing approaches using short-read data can improve consensus accuracy, alignment ambiguity in repetitive regions often leaves some errors unresolved, underscoring the need to balance assembly contiguity, accuracy, and computational complexity when selecting an assembly strategy [5].

Following assembly, alignment to a reference genome provides a framework for assessing assembly quality and identifying genomic variation. Reference-based comparison enables efficient detection of single nucleotide polymorphisms and small insertions or deletions when a closely related reference genome is available [6]. 

In this analysis, raw ONT R10 sequencing reads from Salmonella enterica will be assembled into a consensus genome and compared to a reference genome obtained from NCBI. Several assembly and analysis tools were evaluated based on their suitability for long-read bacterial genome assembly, computational efficiency, and accuracy. The selected workflow balances contiguity, error correction, and interpretability, and includes genome assembly, quality assessment, reference alignment, variant calling, and visualization to evaluate both assembly performance and genomic variation.

Several long-read assemblers were evaluated for assembling the Salmonella enterica genome from Oxford Nanopore data, including Flye, Canu, and NECAT. Flye was selected as the primary assembler due to its robustness and consistent performance for long-read de novo bacterial genome assembly [7,8]. Flye constructs an assembly graph by connecting error-prone disjoint genomic segments into a unified structure, enabling accurate reconstruction of complex genomic regions, and has been shown to perform particularly well for plasmid assembly [7]. While Canu applies extensive read correction and trimming to reduce base-calling errors and can generate reliable assemblies, it has substantially longer runtimes and higher computational demands, making it less efficient for this analysis [7,8]. NECAT employs a two-stage correction and assembly strategy designed to address Nanopore-specific errors and can generate high-quality assemblies more rapidly than Canu, however it often requires additional polishing steps to achieve comparable accuracy [7]. Based on robustness, accuracy, and computational efficiency, Flye was chosen as the most suitable assembler for this workflow.

For reference alignment, minimap2 was chosen due to its widespread use and optimized performance for mapping long-read assemblies to closely related reference genomes [9]. Variant detection and file processing will be conducted using SAMtools and BCFtools for manipulating alignment files and identifying genomic variation. Finally, the Integrative Genomics Viewer (IGV) will be used for visualization, allowing manual inspection of alignments and variants to assess assembly quality and validate detected sequence differences. Together, these tools form a workflow that balances computational efficiency, accuracy, and interpretability for long-read bacterial genome assembly and comparison.

# Methods
## Computational Resources

All analyses were performed on the Digital Research Alliance of Canada’s Nibi high-performance computing cluster. Data processing and quality control steps were executed in an interactive shell environment using the StdEnv/2023 software stack. Containerized software was executed using Apptainer, while select sequence-processing tools were accessed via environment modules provided by the cluster. No raw sequencing data were stored in the version-controlled repository; only scripts and processed result files were archived.

## Sequencing Data Acquisition

Oxford Nanopore Technologies (ONT) R10.4 sequencing reads for Salmonella enterica were obtained from the NCBI Sequence Read Archive (SRA) under accession SRR32410565. Raw sequencing data were downloaded using the SRA Toolkit and compressed in FASTQ format for downstream analysis.

## Read Quality Assessment and Filtering
### Pre-filtering Quality Assessment

Initial read quality was assessed prior to filtering using NanoPlot v1.46.2, executed within an Apptainer container. NanoPlot was used to evaluate read length distributions, quality score distributions, and sequencing yield. Quality control results were generated as interactive HTML reports to enable visual inspection of raw read characteristics.

### Read Filtering

Reads were filtered to improve assembly accuracy by removing low-quality and short reads. Filtering was performed using SeqKit v2.5.1, accessed via the Nibi module system. Reads shorter than 1,000 base pairs or with a mean Phred quality score below 10 were excluded. The resulting filtered reads were written to a new FASTQ file for downstream analysis.

### Post-filtering Quality Assessment

Following filtering, read quality was reassessed using NanoPlot v1.46.2 with the same parameters as the pre-filtering assessment. Post-filtering quality control reports were generated to allow direct comparison of read length, quality distributions, and sequencing yield before and after filtering.

## Genome Assembly

Filtered ONT reads were assembled de novo using Flye v2.9.6, a long-read assembler optimized for error-prone Nanopore sequencing data. Assembly was performed in high-accuracy Nanopore mode using the --nano-hq option, with an estimated genome size of 4.8 Mb, consistent with the expected genome size of Salmonella enterica. Flye was executed using an Apptainer container, and internal polishing steps implemented by Flye were applied automatically to generate a consensus assembly. Assembly outputs included the final assembled contigs and summary assembly statistics.

## Assembly Output Inspection

The resulting assembly consisted of multiple contigs, consistent with a bacterial chromosome and potential plasmid sequences. Assembly statistics, including contig number and length metrics, were obtained from the Flye assembly summary files for downstream quality assessment and comparison to a reference genome.

## Assembly Quality Assessment
Assembly quality was evaluated using QUAST v5.3.0, which compares assembled contigs to a reference genome to assess contiguity, completeness, and structural accuracy. The Flye-assembled genome was compared to the Salmonella enterica reference genome ASM694v2 obtained from NCBI. QUAST was used to calculate assembly metrics including total assembly length, number of contigs, N50, GC content, genome fraction, and the number of misassemblies relative to the reference genome. These metrics were used to assess overall assembly quality prior to downstream reference alignment and variant analysis.

## Assembly Polishing and Post-polishing Quality Assessment

The Flye assembly was further polished using Medaka, which performs neural-network-based consensus correction optimized for Oxford Nanopore sequencing data. Polishing was conducted using filtered ONT R10.4 reads aligned to the Flye assembly, with the Medaka model r1041_e82_400bps_sup_v5.2.0, selected for high-accuracy Nanopore reads. Medaka was executed using an Apptainer container to generate a polished consensus assembly.

To assess the impact of polishing on assembly quality, the Medaka-polished assembly was re-evaluated using QUAST v5.3.0 by comparison to the Salmonella enterica reference genome (ASM694v2). Assembly metrics including total length, contig count, N50, GC content, genome fraction, and misassemblies were compared to the pre-polishing QUAST results to evaluate changes resulting from consensus correction.

## Assembly-to-Reference Alignment

To evaluate the structural concordance between the assembled genome and the reference sequence, the Medaka-polished assembly was aligned to the Salmonella enterica reference genome (ASM694v2) using Minimap2 with the asm5 preset, which is optimized for assembly-to-reference alignment assuming moderate sequence divergence. The resulting SAM file was converted to BAM format, sorted, and indexed using SAMtools, enabling visualization and assessment of contig placement relative to the reference genome

## Read Alignment to Reference Genome

Filtered Oxford Nanopore R10.4 reads were aligned to the Salmonella enterica reference genome (ASM694v2) using Minimap2 with the map-ont preset, which is optimized for long-read Oxford Nanopore sequencing data. The resulting alignments were converted to BAM format, sorted, and indexed using SAMtools, generating alignment files suitable for downstream variant calling and visualization.

## Variant Calling

Variant calling was performed using Clair3, a neural-network-based variant caller optimized for Oxford Nanopore sequencing data. Prior to variant calling, the reference genome was indexed using SAMtools to generate a FASTA index required by Clair3. The sorted and indexed read alignments were provided as input, and variants including single nucleotide polymorphisms and small insertions and deletions were identified relative to the reference genome using a model trained on high-accuracy ONT R10.4 reads.

# Genome assembly and polishing

Oxford Nanopore R10 sequencing reads (FASTQ format) will be quality-checked using NanoPlot (v1.46.2) to assess read length distributions, quality scores, and sequencing yield, ensuring sufficient coverage and read length for reliable _de novo_ assembly.
Following the QC steps, reads will be assembled using Flye (v2.9.6), a long-read assembler designed for error-prone Nanopore data. Assembly was performed in high-accuracy Nanopore mode `--nano-hq`, with an estimated genome size of 4.8 Mb `--genome-size 4.8m`, consistent with _Salmonella enterica_ [10]. Default parameters will be used unless otherwise specified. Flye’s internal polishing steps will be applied to improve consensus accuracy prior to downstream analysis.

`flye --nano-hq reads.fastq --genome-size 4.8m --out-dir flye_output --threads 8`

## Reference genome alignment
To evaluate the assembled genome and identify sequence differences, the draft assembly will be aligned to a _Salmonella enterica_ reference genome downloaded from the NCBI RefSeq database. Alignment will be performed using minimap2 (v2.26), using the `-ax asm5` option to align assembled contigs to the reference, output to a SAM format, assuming ~5% divergence.

`minimap2 -ax asm5 reference.fasta flye_output/assembly.fasta > alignment.sam`

## Alignment processing and variant calling
Alignment files (SAM files) will be converted to BAM format, sorted, and indexed using SAMtools (v1.19.2) to prepare them for variant analysis. Variants relative to the reference genome, including single nucleotide polymorphisms (SNPs) and small insertions and deletions (indels), will be identified using BCFtools (v1.19).

## Visualization
To visually inspect alignment quality and genomic differences between the assembled genome and the reference, alignments and variant calls will be examined using Integrative Genomics Viewer (IGV, v2.19.7). 

# Results

|Metric|Raw Reads|Filtered Reads|
|------|---------|--------------|
|Mean read length|4,128.4|4,165.7|
|Mean read quality|18.9|20.3|
|Median read length|3,957.0|3,980.0|
|Median read quality|23.7|23.8|
|Number of reads|196,031.0|188,542.0|
|Read length N50|4,683.0|4,688.0|
|Total bases|809,296,219.0|785,406,207.0|

Table 1. Summary statistics of Oxford Nanopore R10.4 sequencing reads before and after filtering, showing read length, quality, and yield metrics for the raw dataset and for reads retained after filtering using SeqKit to remove sequences shorter than 1 kb and with average Phred quality scores below 10.

<p align="center">
  <img src="https://github.com/user-attachments/assets/5608a824-9667-40d8-81bf-4ae398457017" width="45%" />
  <img src="https://github.com/user-attachments/assets/b32452a6-4fd1-4ac9-810b-d1aaad8d8eaa" width="45%" />
</p>

Figure 1. NanoPlot-generated KDE plots illustrating read length and average read quality distributions before and after filtering, demonstrating improved read quality and retention of long-read characteristics following filtering.

Quality control analysis of the raw Oxford Nanopore R10.4 reads showed a broad distribution of read lengths centered around ~4 kb, with a median read quality of 23.7 and an N50 of 4,683 bp (Table 1). While the majority of reads exhibited high quality, the KDE plot of read length versus average read quality revealed a distinct low-quality tail, particularly among shorter reads (Figure 1). This skew is reflected by the mean read quality (18.9) being substantially lower than the median, indicating that a subset of lower-quality reads reduced the overall average. Following filtering with SeqKit to remove reads shorter than 1 kb and with average Phred quality scores below 10, the overall read length distribution was largely preserved, as evidenced by similar mean, median, and N50 read lengths (Table 1). In contrast, average read quality increased after filtering, and the KDE distribution shifted toward higher-quality reads, with reduced density in the low-quality region (Figure 1). These results demonstrate that filtering effectively removed lower-quality reads while retaining long-read characteristics and sufficient sequencing depth for downstream genome assembly.

<img width="828" height="88" alt="image" src="https://github.com/user-attachments/assets/8d13c962-839c-4253-b90d-78e44bcf83d8" />

Table 2. Summary of contigs produced by de novo assembly with Flye. The assembly consists of three contigs totaling approximately 5.1 Mb, including two large contigs representing the bacterial chromosome and one smaller circular contig consistent with a plasmid, indicating a highly contiguous assembly typical of bacterial genomes.

The de novo assembly produced three contigs with lengths of approximately 3.3 Mb, 1.7 Mb, and 0.1 Mb, together accounting for the expected genome size of Salmonella enterica. The presence of only a small number of contigs indicates a high level of assembly contiguity, suggesting that long Oxford Nanopore reads effectively spanned repetitive regions that commonly fragment short-read assemblies. The two large contigs likely represent the chromosomal backbone, while the smaller circular contig is consistent with a plasmid, a common feature of Salmonella genomes. The large contig sizes and high sequencing coverage support the reliability of the assembly and indicate that the genome structure was largely reconstructed without extensive fragmentation.

<img width="2048" height="1382" alt="image" src="https://github.com/user-attachments/assets/fe88139a-a2e0-4548-aea1-ec07d0a6f91e" />
Figure 3. Cumulative contig length comparison between the Medaka-polished assembly and the reference genome generated by QUAST.
Contigs are ordered from largest to smallest, and cumulative assembly length is shown for the polished consensus assembly (red) relative to the reference genome (dashed). The assembly reaches approximately 5.1 Mb across three contigs, closely matching the reference genome length and indicating high contiguity, with most of the genome represented by the two largest contigs.

|Metric|Assembled|Polished|
|------|---------|--------------|
|Genome fraction (%)|95.669|95.669|
|Largest contig|3318776|3318770|
|N50|3318776|3318770|
|GC (%)|52.19|52.19|
|# contigs|3|3|
|# mismatches per 100 kbp|27.37|27.1|
|# indels per 100 kbp|3.79|3.62|
|# misassemblies|25 (relocations)|25 (relocations)|

Table 2. QUAST assembly quality metrics before and after Medaka polishing.
The table summarizes reference-based assembly statistics for the Flye-assembled genome prior to polishing (“Assembled”) and after Medaka polishing (“Polished”), including genome completeness, contiguity, GC content, and structural accuracy relative to the Salmonella enterica reference genome.

Assembly quality metrics before and after Medaka polishing indicate that the overall structure and completeness of the genome remained stable, with an identical genome fraction (95.7%), three contigs, and nearly unchanged largest contig size and N50 values (Table X). In contrast, base-level accuracy metrics showed modest improvement following polishing, with reductions in both mismatches and indels per 100 kb, reflecting correction of residual sequencing errors common in Oxford Nanopore data. The low and consistent GC content further supports assembly stability, while the unchanged number of misassemblies suggests that polishing primarily improved consensus accuracy rather than altering large-scale genome structure. Together, these metrics indicate that Medaka polishing refined nucleotide-level accuracy while preserving the high contiguity of the original Flye assembly.

<p align="center">
  <img width="2998" height="1896" alt="image" src="https://github.com/user-attachments/assets/8c67d668-38e1-4298-804c-04e981081c8b" />
  <img width="2996" height="1894" alt="image" src="https://github.com/user-attachments/assets/5cefff82-51ec-4131-a729-b7751c06846f" />
</p>

Figure 4. IGV visualizations of representative single-nucleotide variants identified by Clair3.
(A) IGV visualization of a synonymous single-nucleotide substitution (T → C) at position NC_003197.2:1,165,215 in the reference genome. The variant is supported by multiple aligned reads and does not alter the encoded amino acid, with the codon continuing to encode valine (V), indicating a synonymous mutation.
(B) IGV visualization of a nonsynonymous single-nucleotide substitution (A → G) at position NC_003197.2:3,113,398. This variant results in an amino acid change from isoleucine (I) to valine (V), and is supported by consistent read evidence across the alignment, indicating a true biological variant rather than a sequencing artifact
Figure 3: T -> C (transition mutation (point)), amino acid remains to be valine




# References
[1]  Baker, M. (2012). De novo genome assembly: what every biologist should know. Nature Methods, 9(4), 333–337. https://doi.org/10.1038/nmeth.1935

[2] Zhou, X., & Faust, K. (2025). A high-throughput and time-efficient Nanopore full-length 16S rRNA gene sequencing protocol for synthetic microbial communities. Methods (San Diego, Calif.), 240, 14–20. https://doi.org/10.1016/j.ymeth.2025.04.003

[3] Altermann, E., Tegetmeyer, H. E., & Chanyi, R. M. (2022). The evolution of bacterial genome assemblies—Where do we need to go next. Microbiome Research Reports, 1(2), 15. https://doi.org/10.20517/mrr.2022.02

[4] Bogaerts, B., Maex, M., Commans, F., Goeders, N., Van den Bossche, A., De Keersmaecker, S. C. J., Roosens, N. H. C., Ceyssens, P.-J., Mattheus, W., & Vanneste, K. (2025). Oxford Nanopore Technologies R10 sequencing enables accurate cgMLST-based bacterial outbreak investigation of Neisseria meningitidis and Salmonella enterica when accounting for methylation-related errors. Journal of Clinical Microbiology, 63(10), e0041025. https://doi.org/10.1128/jcm.00410-25

[5] Wick, R. R., & Holt, K. E. (2022). Polypolish: Short-read polishing of long-read bacterial genome assemblies. PLoS Computational Biology, 18(1), e1009802. https://doi.org/10.1371/journal.pcbi.1009802

[6] Vo, N. S., Tran, Q., Niraula, N., & Phan, V. (2014). RandAL: a randomized approach to aligning DNA sequences to reference genomes. BMC Genomics, 15(Suppl 5), Article S2. https://doi.org/10.1186/1471-2164-15-S5-S2

[7] Safar, H. A., Alatar, F., Nasser, K., Al-Ajmi, R., Alfouzan, W., & Mustafa, A. S. (2023). The impact of applying various de novo assembly and correction tools on the identification of genome characterization, drug resistance, and virulence factors of clinical isolates using ONT sequencing. BMC Biotechnology, 23(1), Article 26. https://doi.org/10.1186/s12896-023-00797-3

[8] Wick, R. R., & Holt, K. E. (2019). Benchmarking of long-read assemblers for prokaryote whole genome sequencing [version 1; peer review: 4 approved]. F1000 Research, 8, 2138. https://doi.org/10.12688/f1000research.21782.1

[9] Li, H. (2021). New strategies to improve minimap2 alignment accuracy. Bioinformatics, 37(23), 4572–4574. https://doi.org/10.1093/bioinformatics/btab705

[10] U.S. National Library of Medicine. (n.d.). Salmonella enterica genome assembly ASM250787v2 - NCBI - NLM. National Center for Biotechnology Information. https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002507875.2/ 

[11] Mikolmogorov. (n.d.). Flye/Docs/USAGE.md at flye · Mikolmogorov/Flye. GitHub. https://github.com/mikolmogorov/Flye/blob/flye/docs/USAGE.md 

[12] Li, H. (n.d.). LH3/Minimap2: A versatile pairwise aligner for genomic and spliced nucleotide sequences. GitHub. https://github.com/lh3/minimap2?tab=readme-ov-file#general 
