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

---

### 1. Sequencing data acquisition and characteristics
Oxford Nanopore long-read sequencing data for *Salmonella enterica* were obtained from the NCBI Sequence Read Archive (SRA) under accession **SRR32410565** :contentReference[oaicite:1]{index=1}. Reads were provided in FASTQ format and generated using R10.4.1 chemistry.

Because ONT datasets typically contain heterogeneous read lengths and variable per-read quality, preprocessing steps were performed prior to genome assembly to remove short or low-confidence reads that could negatively influence graph construction and polishing.

---

### 2. Read quality control and filtering

### Initial quality assessment
Raw ONT reads were first evaluated using NanoPlot to characterize read length distributions, sequencing yield, and per-read quality metrics.

```
NanoPlot --fastq SRR32410565.fastq -o nanoplot_output
```

The generated HTML report (NanoPlot-report.html) was inspected to determine whether the dataset contained sufficient long, high-quality reads to support de novo assembly.

### Read filtering

Reads were then filtered using Filtlong to preferentially retain longer and higher-quality sequences while discarding reads likely to introduce assembly errors. The following thresholds were applied:

- minimum read length: 1000 bp

- retain the best 90% of reads by quality

- minimum mean read quality: Q20

```
filtlong --min_length 1000 --keep_percent 90 -q 20 \
  SRR32410565.fastq > SRR32410565_q20.fastq
```
Filtering reduced low-quality tails and enriched for longer reads suitable for accurate repeat resolution during assembly.

---

### 3. Genome assembly
Filtered reads were assembled de novo using Flye, a repeat-graph–based assembler optimized for long, error-prone nanopore sequencing data.

```
flye -t 4 \
     --genome-size 4.5m \
     --asm-coverage 100 \
     --nano-hq SRR32410565_q20.fastq \
     -o flye_output_q20
```

```--genome-size 4.5m``` provides an approximate genome length for Salmonella enterica, enabling Flye to optimize graph construction and coverage expectations.

```--asm-coverage 100``` instructs Flye to use the longest reads up to ~100× coverage for initial repeat graph construction, improving computational efficiency while preserving assembly quality.

```--nano-hq``` is appropriate for high-accuracy ONT chemistries and enables internal error correction tuned for improved base quality.

Assembly metrics and contig information were obtained from Flye output files, including statistics describing total assembly length, number of contigs, N50, and coverage estimates, which are used for downstream evaluation and comparison.


### 3.1 Assembly polishing
Although long-read assemblers generate highly contiguous genomes, residual base-level errors may remain due to systematic sequencing inaccuracies. To improve consensus accuracy, the Flye assembly was polished using Medaka, a neural network–based polishing tool developed for Oxford Nanopore data.

```
medaka_consensus \
  -i SRR32410565.fastq \
  -d flye_output_q20/assembly.fasta \
  -o medaka_polish \
  -t 4 \
  -m r1041_e82_400bps_sup_v5.0.0
```
The resulting polished consensus sequence (medaka_polish/consensus.fasta) was used for subsequent reference comparisons, quality assessment, and variant analysis.

---

### 4. Reference-based alignment and file processing

### Reference genome
The *Salmonella enterica* reference genome (accession GCF_000006945.2) was downloaded from NCBI using the datasets command-line tool.

### Read-to-reference alignment

To evaluate mapping performance, coverage uniformity, and enable variant calling, filtered reads were aligned to the reference genome using minimap2 with parameters optimized for ONT data.
```
minimap2 -ax map-ont -t 4 \
  ncbi_dataset/data/GCF_000006945.2/GCF_000006945.2_ASM694v2_genomic.fna \
  SRR32410565.fastq > reads_vs_ref.sam
```

SAM output was converted into compressed, sorted, and indexed BAM format to enable efficient querying and visualization.

```
samtools view -bS reads_vs_ref.sam > reads_vs_ref.bam
samtools sort reads_vs_ref.bam -o reads_vs_ref.sorted.bam
samtools index reads_vs_ref.sorted.bam
```

### Assembly-to-reference alignment

To assess structural concordance between the polished assembly and the reference genome, the Medaka consensus was aligned using minimap2 with assembly-aware presets.
```
minimap2 -ax asm5 -t 4 \
  ncbi_dataset/data/GCF_000006945.2/GCF_000006945.2_ASM694v2_genomic.fna \
  medaka_polish/consensus.fasta > assembly_vs_ref.sam
```

SAM output was similarly converted into compressed, sorted, and indexed BAM format to enable efficient querying and visualization.

```
samtools view -bS assembly_vs_ref.sam > assembly_vs_ref.bam
samtools sort assembly_vs_ref.bam -o assembly_vs_ref.sorted.bam
samtools index assembly_vs_ref.sorted.bam
```

---

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





