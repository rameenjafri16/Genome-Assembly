# Genome-Assembly

## Table of Contents
- [General Overview](#general-overview)
- [Introduction](#introduction)
- [Proposed Methods](#proposed-methods)
  - [Sequencing data acquisition and characteristics](#1-sequencing-data-acquisition-and-characteristics)
  - [Read quality control and filtering](#2-read-quality-control-and-filtering)
  - [Genome assembly](#3-genome-assembly)
  - [Reference-based alignment and file processing](#4-reference-based-alignment-and-file-processing)
  - [Visualization](#5-visualization)
- [Citations](#citations)

## General Overview 
This project processes long-read sequencing data generated using Oxford Nanopore Technologies (ONT) to assemble the genome of Salmonella enterica. Raw sequencing reads are first inspected and filtered to assess read length and quality, ensuring that low-quality and short reads are removed prior to assembly. The filtered reads are then assembled de novo using the Flye assembler, producing a highly contiguous genome assembly.

To evaluate assembly quality, the resulting contig is compared to a reference genome obtained from the National Center for Biotechnology Information (NCBI) using minimap2 for long-read alignment. Alignment outputs are processed and visualized using genome browser tools to assess coverage patterns and mapping consistency across the assembly. Together, this workflow integrates quality control, assembly, alignment, and visualization to evaluate the feasibility, advantages, and limitations of long-read genome assembly and reference-based comparison for bacterial genomes.


## Introduction
Long-read sequencing became widely accessible with the release of the MinION by Oxford Nanopore Technologies (ONT) in 2014, enabling real-time sequencing of native DNA and RNA (Peng et al., 2025). Since then, additional long-read platforms have emerged, including Pacific Biosciences (PacBio) and high-throughput systems such as QitanTech and CycloneSEQ (Peng et al., 2025). These platforms differ in accuracy, cost, and accessibility, influencing their suitability for routine bacterial genome assembly (Peng et al., 2025). PacBio HiFi sequencing can achieve base-level accuracies exceeding 99.99% but requires expensive instrumentation and infrastructure, limiting its use in smaller laboratories and field-based applications (Peng et al., 2025). Other long-read platforms can generate reads of sufficient length for genome assembly but often exhibit lower base-level accuracy, with reported average read lengths of approximately 11.6 kb and quality scores around Q14.4 (Peng et al., 2025).

Oxford Nanopore sequencing provides a balance of read length, accessibility, and cost that makes it well suited for bacterial genome assembly (Wang et al., 2021). Unlike short-read sequencing platforms such as Illumina or MGI, which produce reads of only a few hundred base pairs, nanopore sequencing generates long reads capable of spanning multidrug resistance (MDR) regions and other complex genomic structures (Peng et al., 2025). This is particularly important for bacterial genomes, where antimicrobial resistance genes (ARGs) and mobile genetic elements such as plasmids, transposons, and integrons often occur within repetitive or structurally complex regions (Peng et al., 2025). In addition, ONT sequencing enables real-time data generation, portability, and relatively low sequencing costs, while barcoding strategies allow multiple bacterial genomes to be sequenced simultaneously on a single flow cell, reducing per-genome sequencing costs and supporting rapid antimicrobial resistance surveillance (Peng et al., 2025).

Early ONT sequencing technologies were limited by raw read accuracies of approximately 90%, constraining their use for accurate genome assembly (Peng et al., 2025). Continued improvements in nanopore chemistry and basecalling algorithms have substantially increased sequencing accuracy (Wang et al., 2024). The introduction of the R10 nanopore, featuring a longer barrel and dual reader head, improved basecalling accuracy in homopolymeric regions, while subsequent R10.4 and R10.4.1 flow cells and Q20+ chemistry enabled raw read accuracies exceeding 99% (Zhang et al., 2023). These advances have made long-read-only microbial genome assembly increasingly feasible (Zhang et al., 2023).

Despite improved platform accuracy, long-read datasets remain heterogeneous in read length and per-read quality, which can negatively impact genome assembly if not addressed (Pardo-Palacios, 2024). Low-quality or short reads may introduce assembly errors and reduce contiguity, making read-level quality assessment and filtering critical preprocessing steps in long-read genome assembly workflows (Jiao et al., 2017).

Accurate comparison of an assembled genome to a reference requires effective alignment of long, error-prone reads or contigs (Jiao et al., 2017). Alignment tools originally designed for Sanger or short-read sequencing data are poorly suited for long-read sequencing due to higher error rates and read lengths (Kent, 2002; Wang et al., 2021). Specialized long-read aligners were therefore developed, and as ONT read lengths increased beyond 100 kb, minimap2 was introduced using a seed–chain–align strategy optimized for long-read data (Wang et al., 2021). Benchmarking studies have shown that minimap2 achieves faster runtimes than alternative long-read aligners without sacrificing accuracy, making it well suited for reference-based alignment of genomic ONT data (Wang et al., 2021). Visualization of alignment results provides an additional qualitative evaluation step, allowing assessment of coverage patterns, mapping consistency, and potential assembly artifacts.

Overall, advances in ONT sequencing chemistry, read preprocessing, and long-read alignment algorithms have made long-read-only bacterial genome assembly increasingly viable (Wang et al., 2021). However, outcomes remain sensitive to parameter choices, residual sequencing errors, and reference genome quality. In this study, these approaches are applied to assemble and align the genome of Salmonella enterica to evaluate the strengths and limitations of long-read genome assembly and reference-based comparison (Wang et al., 2021).



## Proposed Methods
### Sequencing data acquisition and characteristics
Oxford Nanopore long-read sequencing data for *Salmonella enterica* were obtained from the NCBI Sequence Read Archive (SRA) under accession **SRR32410565** :contentReference[oaicite:1]{index=1}. Reads were provided in FASTQ format and generated using R10.4.1 chemistry.

Because ONT datasets typically contain heterogeneous read lengths and variable per-read quality, preprocessing steps were performed prior to genome assembly to remove short or low-confidence reads that could negatively influence graph construction and polishing.


### 2. Read quality control and filtering
Read quality was assessed both before and after filtering using NanoPlot, which generates summary statistics and an HTML report describing read length distributions and per-read quality patterns for long-read sequencing data. NanoPlot was used descriptively to evaluate whether the raw dataset contained sufficient read length and quality for de novo genome assembly, and to confirm that filtering improved overall read quality.

Reads were filtered using Filtlong to remove short and low-quality reads. Filtering thresholds were chosen to retain informative long reads while discarding reads likely to introduce errors during assembly. Reads shorter than 1000 bp were removed (--min_length 1000), only the top 90% of reads ranked by quality were retained (--keep_percent 90), and reads with mean quality scores below Q20 were excluded (-q 20). The filtered FASTQ output was then re-evaluated with NanoPlot to confirm improved read quality and a clear cutoff corresponding to the applied thresholds.

### 3. Genome assembly
Filtered reads were assembled de novo using Flye with the --nano-hq preset, which is intended for higher-accuracy nanopore reads. Flye was selected because it is optimized for long-read assembly and can generate contiguous assemblies from error-prone long reads. The assembly produced a single contig, consistent with a highly contiguous bacterial genome assembly suitable for downstream comparison to a reference genome.

### 4. Reference-based alignment and file processing
To compare the assembly to sequencing reads and support downstream inspection, alignments were generated using minimap2. Minimap2 was run using the ONT mapping preset and configured to output alignments in SAM format using 32 threads. The resulting SAM file was intended for conversion into compressed, indexed formats for visualization. SAM outputs were converted to BAM, sorted, and indexed using samtools to ensure compatibility with genome browsers and efficient navigation across the genome.

### 5. Visualization
Alignments were visualized using the Integrative Genomics Viewer (IGV). The assembled genome and the sorted, indexed BAM file were loaded into IGV to inspect alignment consistency and coverage patterns across the contig. Visualization enabled qualitative assessment of whether reads mapped cleanly to the assembly and whether any regions showed unusual coverage patterns that could indicate assembly artifacts or difficult-to-map regions.



## Citations 
Zhang, T., Li, H., Ma, S., Cao, J., Liao, H., Huang, Q., & Chen, W. (2023). The newest Oxford Nanopore R10.4.1 full-length 16S rRNA sequencing enables the accurate resolution of species-level microbial community profiling. Applied and environmental microbiology, 89(10), e0060523. https://doi.org/10.1128/aem.00605-23

Oxford Nanopore Technologies. (2025, February 19). How Oxford Nanopore sequencing works. Oxford Nanopore Technologies. Retrieved 2026, January 16, from https://nanoporetech.com/blog/how-oxford-nanopore-sequencing-works

Peng, K., Li, C., Wang, Q. et al. The applications and advantages of nanopore sequencing in bacterial antimicrobial resistance surveillance and research. npj Antimicrob Resist 3, 87 (2025). https://doi.org/10.1038/s44259-025-00157-5 

Wang, Y., Zhao, Y., Bollas, A. et al. Nanopore sequencing technology, bioinformatics and applications. Nat Biotechnol 39, 1348–1365 (2021). https://doi.org/10.1038/s41587-021-01108-x 

Pardo-Palacios, F.J., Wang, D., Reese, F. et al. Systematic assessment of long-read RNA-seq methods for transcript identification and quantification. Nat Methods 21, 1349–1363 (2024). https://doi.org/10.1038/s41592-024-02298-3

Wang, Z., Fang, Y., Liu, Z. et al. Adapting nanopore sequencing basecalling models for modification detection via incremental learning and anomaly detection. Nat Commun 15, 7148 (2024). https://doi.org/10.1038/s41467-024-51639-5

Jiao, W. B., Accinelli, G. G., Hartwig, et al. Improving and correcting the contiguity of long-read genome assemblies of three plant species using optical mapping and chromosome conformation capture data. Genome research (2017), 27(5), 778–786. https://doi.org/10.1101/gr.213652.116 

Kent, W James. “BLAT--the BLAST-like alignment tool.” Genome research vol. 12,4 (2002): 656-64. doi:10.1101/gr.229202



## Software & Resources
**Core tools**
- Flye (assembly): https://github.com/fenderglass/Flye
- minimap2 (alignment): https://github.com/lh3/minimap2
- samtools (SAM/BAM processing): https://github.com/samtools/samtools

**QC + filtering**
- NanoPlot (QC reporting): https://github.com/wdecoster/NanoPlot
- Filtlong (read filtering): https://github.com/rrwick/Filtlong

**Visualization**
- IGV (genome browser): https://software.broadinstitute.org/software/igv/

**Data source**
- NCBI SRA record (SRR32410565): https://www.ncbi.nlm.nih.gov/sra/SRR32410565

## Software Versions & Reproducibility
All software was installed via **bioconda** and executed within a conda environment on an **Ubuntu virtual machine**.
All analyses were performed using open-source software. Software versions and parameters are documented to support reproducibility of the workflow.

| Tool | Version |
|-----|--------|
| Flye | 2.9.6-b1802 |
| NanoPlot | 1.46.1 |
| Filtlong | 0.3.1 |
| minimap2 | 2.30-r1287 |
| samtools | 1.22.1 |





