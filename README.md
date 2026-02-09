# Genome-Assembly

## Table of Contents
- [General Overview](#general-overview)
- [Introduction](#introduction)
- [Proposed Methods](#proposed-methods)
  - [Sequencing data acquisition and characteristics](#1-sequencing-data-acquisition-and-characteristics)
  - [Read quality control and filtering](#2-read-quality-control-and-filtering)
  - [Genome assembly](#3-genome-assembly)
  - [Assembly polishing](#31-assembly-polishing)
  - [Reference-based alignment and file processing](#4-reference-based-alignment-and-file-processing)
  - [Variant calling](#5-variant-calling)
  - [Visualization](#6-visualization)
  - [Coverage Analysis](#coverage-analysis)
  - [Variant Density Visualization](#variant-density-visualization)
- [Results](#results)
  - [Assembly Quality and Alignment to Reference](#assembly-quality-and-alignment-to-reference)
  - [Read Alignment and Variant Analysis](#read-alignment-and-variant-analysis)
- [Discussion](#discussion)
  - [Assembly Quality and Workflow Success](#assembly-quality-and-workflow-success)
  - [Plasmid Divergence and Biological Implications](#plasmid-divergence-and-biological-implications)
- [Citations](#citations)
- [Software & Resources](#software--resources)
- [Software Versions & Reproducibility](#software-versions--reproducibility)


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



## Methods

---

### 1. Sequencing data acquisition and characteristics
Oxford Nanopore long-read sequencing data for *Salmonella enterica* were obtained from the NCBI Sequence Read Archive (SRA) under accession **SRR32410565** :contentReference[oaicite:1]{index=1}. Reads were provided in FASTQ format and generated using R10.4.1 chemistry.

Because ONT datasets typically contain heterogeneous read lengths and variable per-read quality, preprocessing steps were performed prior to genome assembly to remove short or low-confidence reads that could negatively influence graph construction and polishing.

---

### 2. Read quality control and filtering

#### Initial quality assessment
Raw ONT reads were first evaluated using NanoPlot to characterize read length distributions, sequencing yield, and per-read quality metrics.

```
NanoPlot --fastq SRR32410565.fastq -o nanoplot_output
```

The generated HTML report (NanoPlot-report.html) was inspected to determine whether the dataset contained sufficient long, high-quality reads to support de novo assembly.

#### Read filtering

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

#### Reference genome
The *Salmonella enterica* reference genome (accession GCF_000006945.2) was downloaded from NCBI using the datasets command-line tool.

#### Read-to-reference alignment

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

#### Assembly-to-reference alignment

To assess structural concordance between the polished assembly and the reference genome, the Medaka consensus was aligned using minimap2 with assembly-aware presets.
```
minimap2 -ax asm5 -t 4 \
  ncbi_dataset/data/GCF_000006945.2/GCF_000006945.2_ASM694v2_genomic.fna \
  medaka_polish/consensus.fasta > assembly_vs_ref.sam
```

SAM output was similarly converted into compressed, sorted, and indexed BAM format to enable efficient querying and visualization.

---

### 5. Variant calling
Small variants were identified from the read-to-reference alignments using Clair3, a deep learning–based variant caller optimized for long-read sequencing data.

#### Reference preparation
Prior to variant calling, the reference FASTA was indexed:

```
samtools faidx ncbi_dataset/data/GCF_000006945.2/GCF_000006945.2_ASM694v2_genomic.fna
```

#### Running Clair3

Variant calling was performed using the sorted BAM file generated from minimap2. Analyses were configured for haploid bacterial genomes and executed using ONT-specific models.
``` --platform ont ```  applies ONT-specific error models
``` --include_all_ctgs```  evaluate every reference contig
``` --haploid_precise```  use assumptions appropriate for haploid genomes
``` --no_phasing_for_fa```  disable haplotype phasing

### 6. Visualization
Alignments were visualized using the Integrative Genomics Viewer (IGV). The assembled genome and the sorted, indexed BAM file were loaded into IGV to inspect alignment consistency and coverage patterns across the contig. Visualization enabled qualitative assessment of whether reads mapped cleanly to the assembly and whether any regions showed unusual coverage patterns that could indicate assembly artifacts or difficult-to-map regions.

#### Coverage Analysis
Read coverage across the reference genome was visualized using R (version 4.5.2) with the following packages:

ggplot2 - for creating publication-quality plots
dplyr - for data manipulation
tidyr - for data reshaping
gridExtra - for multi-panel figure layout
scales - for axis formatting

Coverage metrics were extracted from the sorted BAM file using samtools coverage, which provided per-contig statistics including mean depth, coverage percentage, and mapping quality. The R scripts (coverage_simple.R and coverage_analysis.R) generated bar plots comparing coverage metrics between the chromosome (NC_003197.2) and plasmid (NC_003277.2).
*Complete R scripts are available in the repository.*

#### Variant Density Visualization
Genome-wide variant distribution was visualized using a circular plot created with R package circlize and tidyverse. Variants were filtered for quality scores Q≥20 from the Clair3 VCF output. The circular plot displays:
*Script available in the repository.*

## Results 

#### Assembly Quality and Alignment to Reference
The Medaka-polished assembly produced three contigs with a total length of 5,104,809 bp, compared to the reference genome length of 4,951,383 bp. The assembly achieved an N50 of 3,318,771 bp with the largest contig spanning 3,318,771 bp. QUAST analysis showed that 95.67% of the reference genome was covered by the assembly, with a duplication ratio of 1.002. The assembly contained 25 misassemblies across 2 contigs, with 10 local misassemblies and a mismatch rate of 27.10 per 100 kbp. The assembly's GC content was 52.19%, compared to the reference GC content of 52.24%. The largest aligned segment was 953,688 bp, with an NA50 of 460,923 bp.
When the assembled contigs were aligned to the reference genome using minimap2, 96% of assembly contigs mapped successfully. The chromosome (NC_003197.2) showed 97.42% coverage at a mean depth of 0.98x and a mean mapping quality of 60. In contrast, the plasmid (NC_003277.2) showed 0% coverage in the assembly-to-reference alignment, indicating that the assembled plasmid contig does not align to the reference virulence plasmid pSLT.

INSERT 6 PANEL
(Figure 2: Assembly quality metrics and comparison to reference)

#### Read Alignment to Reference Genome
Raw read alignment to the reference genome showed 183,082 reads aligning to the chromosome (NC_003197.2), achieving 97.83% coverage at a mean depth of 150.97x with a mean mapping quality of 59.5. Mean base quality for chromosome-aligned reads was 41.3. The plasmid (NC_003277.2) had 4,186 aligned reads with 45.74% coverage at a mean depth of 98.02x and a mean mapping quality of 51.7. Mean base quality for plasmid-aligned reads was 41.6.
The chromosome achieved nearly complete breadth of coverage (97.83%) at high mean depth (~151×) and mapping quality near maximum, indicating reliable representation of chromosomal sequence in the read data. In contrast, the plasmid showed substantially lower breadth of coverage (45.74%), despite similar read depth in regions where alignment occurred. This pattern suggests that only portions of the plasmid are present or sufficiently similar to the reference, consistent with structural divergence or partial absence.

#### Visual Inspection of Read Alignments
Alignments were inspected using the Integrative Genomics Viewer (IGV) to assess mapping quality and identify potential assembly artifacts. Visual inspection confirmed that reads mapped consistently across most of the chromosome, with few abrupt coverage changes or large gaps. Mismatch patterns appeared dispersed rather than clustered, suggesting that the majority of called variants reflect genuine sequence divergence rather than systematic mapping errors.

For the plasmid, IGV visualization revealed multiple SNPs and short indels supported by consistent read evidence across aligned regions (Figure 1). However, coverage was fragmented with extensive gaps, confirming the quantitative finding of 45.74% reference coverage. Where reads did align, variant density was extremely high, with numerous mismatches visible across short genomic windows.
<img width="1710" height="1107" alt="Screenshot 2026-02-08 at 2 00 22 PM" src="https://github.com/user-attachments/assets/094ef974-94b0-4030-8509-eb152a42df63" />
**Figure 1:** IGV visualization of read alignments to plasmid NC_003277.2. Multiple SNPs and short indels are supported by consistent read evidence across the region, indicating genuine sequence divergence rather than sporadic mapping error. Note the fragmented coverage pattern with gaps between aligned regions.

#### Variant Calling and Genomic Distribution
Variant calling using Clair3 with quality filtering (Q≥20) identified 9,658 variants across the genome, comprising 8,647 SNPs (89.5%), 421 insertions (4.4%), and 590 deletions (6.1%). The distribution of variants differed dramatically between the chromosome and plasmid. The chromosome contained 2,295 variants across 4.86 Mb (0.47 variants per kb), while the plasmid contained 7,363 variants across 0.09 Mb (82.0 variants per kb), representing a 174-fold higher variant density on the plasmid (Figure 4).

Circular visualization of variant density showed relatively uniform SNP distribution across the chromosome with occasional concentration zones (Figure 3). In contrast, plasmid variants clustered heavily in regions that aligned to the reference sequence, with dense variant accumulation visible in 5 kb bins. The extreme variant density in plasmid-aligned regions, combined with the 0% assembly-to-reference alignment and fragmented read coverage, suggests that the assembled plasmid represents a sequence fundamentally different from the reference pSLT virulence plasmid.

<img width="1674" height="918" alt="image" src="https://github.com/user-attachments/assets/21abd60e-f6f6-49b1-9876-343a0082d384" />
**Figure 3:** Circular representation of genome-wide coverage and variant density (Q ≥ 20). Variants are broadly distributed across the chromosome but concentrate heavily within plasmid-associated regions, consistent with elevated divergence and fragmented mapping relative to the reference.

<img width="3418" height="2138" alt="image" src="https://github.com/user-attachments/assets/08d873c6-fdbe-4dbd-96ab-bbd94b9c267b" />
**Figure 4:** Variant density across plasmid NC_003277.2 in 100 bp windows. SNPs dominate the signal and form dense clusters across much of the plasmid, while insertions and deletions occur at lower frequency; coverage is uneven, with several regions showing reduced or absent alignment.


## Discussion 
#### Assembly Quality and Workflow Success
The genome assembly pipeline successfully produced a high-quality draft genome of the Salmonella enterica strain, achieving 95.67% coverage of the reference genome. The N50 of 3.3 Mb and collapse into just three contigs demonstrates that Nanopore long-read sequencing combined with Flye assembly and Medaka polishing effectively resolved the genome structure. Chromosome alignment showed 97.42% coverage with high mapping quality (MAPQ=60), while read-to-reference alignment confirmed adequate sequencing depth (151× mean coverage) and accurate base calling (Q=41.3).
The 2,295 chromosomal variants (0.47 variants/kb) likely represent genuine strain-level differences from the reference strain LT2, consistent with natural genomic variation within Salmonella populations (Robertson et al., 2023). The 25 misassemblies and mismatch rate of 27.10 per 100 kbp are within acceptable ranges for Nanopore-only assemblies and could be further reduced through hybrid assembly with Illumina short reads.

#### Plasmid Divergence and Biological Implications
The most striking finding is the substantial divergence of the assembled plasmid from the reference plasmid pSLT (NC_003277.2). While raw reads showed 45.74% coverage of the reference at 98× depth, the assembled plasmid exhibited 0% alignment to pSLT, indicating a fundamentally different plasmid is present in this strain. This interpretation is strongly supported by the extraordinarily high variant density (82.0 variants/kb)—174-fold higher than the chromosome—suggesting misalignments from forcing reads onto an incorrect reference rather than true polymorphisms.

The reference plasmid pSLT, originally characterized from S. Typhimurium LT2, is a 94-kb conjugative virulence plasmid containing 108 coding sequences, including the spv operon (spvRABCD) that contributes to systemic virulence (McClelland et al., 2001). McClelland et al. demonstrated that pSLT is self-transmissible at rates up to 3 × 10⁻⁴ due to homologues of the F-factor tra operon, and estimated copy numbers of 1.4–3.1 under various growth conditions. Importantly, only three pSLT genes showed homology to other Salmonella serovars (S. Typhi, S. Paratyphi A, S. Paratyphi B), as these strains typically lack this plasmid. However, 50 pSLT genes had close homologues in plasmids from other Salmonella serovars, suggesting a family of related conjugative plasmids circulates within the genus.

The complete absence of alignment between our assembled plasmid and pSLT indicates our strain carries a plasmid from a different incompatibility group or MOB-cluster. This has significant biological implications, as Robertson et al. (2023) identified 1,044 MOB-clusters from 183,017 plasmids across 1,204 Salmonella serovars, with 22% carrying at least one antimicrobial resistance (AMR) gene. Conjugative and mobilizable plasmids show significantly higher serovar diversity than non-mobilizable plasmids (median 2.2 and 2.1 nats vs. 1.8 nats), with 88.3% of broad-host-range plasmids being mobilizable. Experimental work by Laidlaw et al. (2024) confirmed variable conjugation rates (10⁻⁸ to 10⁻⁴) across Salmonella serovars, demonstrating that plasmid transfer dynamics are modulated by both host genetics and environmental factors.

Future characterization should include MOB-cluster typing, replicon typing to determine incompatibility group and host range, functional annotation to identify AMR genes and virulence factors, and conjugation assays to experimentally assess transmissibility. Given ongoing multi-plasmid AMR outbreaks identified by Robertson et al. (2023), understanding whether this alternative plasmid carries clinically relevant determinants could provide valuable insights into horizontal gene transfer dynamics in Salmonella.


## Citations 
Zhang, T., Li, H., Ma, S., Cao, J., Liao, H., Huang, Q., & Chen, W. (2023). The newest Oxford Nanopore R10.4.1 full-length 16S rRNA sequencing enables the accurate resolution of species-level microbial community profiling. Applied and environmental microbiology, 89(10), e0060523. https://doi.org/10.1128/aem.00605-23

Oxford Nanopore Technologies. (2025, February 19). How Oxford Nanopore sequencing works. Oxford Nanopore Technologies. Retrieved 2026, January 16, from https://nanoporetech.com/blog/how-oxford-nanopore-sequencing-works

Peng, K., Li, C., Wang, Q. et al. The applications and advantages of nanopore sequencing in bacterial antimicrobial resistance surveillance and research. npj Antimicrob Resist 3, 87 (2025). https://doi.org/10.1038/s44259-025-00157-5 

Wang, Y., Zhao, Y., Bollas, A. et al. Nanopore sequencing technology, bioinformatics and applications. Nat Biotechnol 39, 1348–1365 (2021). https://doi.org/10.1038/s41587-021-01108-x 

Pardo-Palacios, F.J., Wang, D., Reese, F. et al. Systematic assessment of long-read RNA-seq methods for transcript identification and quantification. Nat Methods 21, 1349–1363 (2024). https://doi.org/10.1038/s41592-024-02298-3

Wang, Z., Fang, Y., Liu, Z. et al. Adapting nanopore sequencing basecalling models for modification detection via incremental learning and anomaly detection. Nat Commun 15, 7148 (2024). https://doi.org/10.1038/s41467-024-51639-5

Jiao, W. B., Accinelli, G. G., Hartwig, et al. Improving and correcting the contiguity of long-read genome assemblies of three plant species using optical mapping and chromosome conformation capture data. Genome research (2017), 27(5), 778–786. https://doi.org/10.1101/gr.213652.116 

Kent, W James. “BLAT--the BLAST-like alignment tool.” Genome research vol. 12,4 (2002): 656-64. doi:10.1101/gr.229202

Laidlaw, A., Blondin-Brosseau, M., Shay, J. A., Dussault, F., Rao, M., Petronella, N., & Tamber, S. (2024). Variation in plasmid conjugation among nontyphoidal Salmonella enterica serovars. Canadian Journal of Microbiology, 70(12), e2024-0164. https://doi.org/10.1139/cjm-2024-0164

McClelland, M., Sanderson, K. E., Spieth, J., Clifton, S. W., Latreille, P., Courtney, L., Porwollik, S., Ali, J., Dante, M., Du, F., Hou, S., Layman, D., Leonard, S., Nguyen, C., Scott, K., Holmes, A., Grewal, N., Mulvaney, E., Ryan, E., Sun, H., Florea, L., Miller, W., Stoneking, T., Nhan, M., Waterston, R., & Wilson, R. K. (2001). Complete genome sequence of Salmonella enterica serovar Typhimurium LT2. Nature, 413(6858), 852–856. https://doi.org/10.1038/35101614

Robertson, J., Schonfeld, J., Bessonov, K., Bastedo, P., & Nash, J. H. E. (2023). A global survey of Salmonella plasmids and their associations with antimicrobial resistance. Microbial Genomics, 9(5), 001002. https://doi.org/10.1099/mgen.0.001002


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


 Tool | Version |
|------|---------|
| Flye | 2.9.6-b1802 |
| Medaka | 1.11.3 |
| NanoPlot | 1.46.1 |
| Filtlong | 0.3.1 |
| minimap2 | 2.30-r1287 |
| samtools | 1.22.1 |
| Clair3 | 1.0.10 |
| QUAST | 5.2.0 |
| IGV | 2.16.2 |
| R | 4.3.1 |
| ggplot2 | 3.4.4 |
| dplyr | 1.1.3 |
| tidyr | 1.3.0 |
| circlize | 0.4.15 |




