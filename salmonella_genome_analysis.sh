#!/bin/bash
# Salmonella enterica Genome Analysis Pipeline
# Nanopore sequencing data analysis workflow

## 01_Download raw reads
# wget https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR32410565&display=download

## 02_Quality assessment of raw reads
nanoplot --fastq SRR32410565.fastq -o nanoplot_output

## 03_Filter reads
filtlong --min_length 1000 --keep_percent 90 -q 20 SRR32410565.fastq > SRR32410565_q20.fastq

## 04_Genome assembly
flye -t 4 --genome-size 4.5m --asm-coverage 100 --nano-hq SRR32410565_q20.fastq -o flye_output_q20

## 05_Assembly polishing
medaka_consensus -i SRR32410565.fastq \
                 -d flye_output_q20/assembly.fasta \
                 -o medaka_polish \
                 -t 4 \
                 -m r1041_e82_400bps_sup_v5.0.0

## 06_Assembly quality assessment
quast.py medaka_polish/consensus.fasta \
         -r ncbi_dataset/data/GCF_000006945.2/GCF_000006945.2_ASM694v2_genomic.fna \
         -o quast_results \
         --threads 4

## 07_Align assembly to reference genome
minimap2 -ax asm5 -t 4 \
  ncbi_dataset/data/GCF_000006945.2/GCF_000006945.2_ASM694v2_genomic.fna \
  medaka_polish/consensus.fasta \
  > assembly_vs_ref.sam

## 08_Align reads to reference genome
minimap2 -ax map-ont -t 4 \
  ncbi_dataset/data/GCF_000006945.2/GCF_000006945.2_ASM694v2_genomic.fna \
  SRR32410565.fastq \
  > reads_vs_ref.sam

## 09_Convert and process reads alignment
samtools view -bS reads_vs_ref.sam > reads_vs_ref.bam
samtools sort reads_vs_ref.bam -o reads_vs_ref.sorted.bam
samtools index reads_vs_ref.sorted.bam
samtools flagstat reads_vs_ref.sorted.bam > reads_alignment_stats.txt
samtools depth reads_vs_ref.sorted.bam > reads_coverage.txt
samtools coverage reads_vs_ref.sorted.bam > reads_coverage_summary.txt

## 10_Convert and process assembly alignment
samtools view -bS assembly_vs_ref.sam > assembly_vs_ref.bam
samtools sort assembly_vs_ref.bam -o assembly_vs_ref.sorted.bam
samtools index assembly_vs_ref.sorted.bam
samtools flagstat assembly_vs_ref.sorted.bam > assembly_alignment_stats.txt
samtools coverage assembly_vs_ref.sorted.bam > assembly_coverage_summary.txt

## 11_Index reference genome for variant calling
samtools faidx ncbi_dataset/data/GCF_000006945.2/GCF_000006945.2_ASM694v2_genomic.fna

## 12_Variant calling with Clair3
docker run --rm --platform linux/amd64 \
  -v "$PWD":"$PWD" -w "$PWD" \
  hkubal/clair3 \
  /opt/bin/run_clair3.sh \
    --bam_fn=reads_vs_ref.sorted.bam \
    --ref_fn=ncbi_dataset/data/GCF_000006945.2/GCF_000006945.2_ASM694v2_genomic.fna \
    --threads=4 \
    --platform=ont \
    --model_path=/opt/models/ont \
    --output=clair3_output \
    --sample_name=Salmonella_sample \
    --include_all_ctgs \
    --haploid_precise \
    --no_phasing_for_fa

## 13_Variant statistics
echo "Total variants:"
gunzip -c clair3_output/merge_output.vcf.gz | grep -v "^#" | wc -l

echo "SNPs:"
gunzip -c clair3_output/merge_output.vcf.gz | grep -v "^#" | \
  awk '{if(length($4)==1 && length($5)==1) print}' | wc -l

echo "Insertions:"
gunzip -c clair3_output/merge_output.vcf.gz | grep -v "^#" | \
  awk '{if(length($4)<length($5)) print}' | wc -l

echo "Deletions:"
gunzip -c clair3_output/merge_output.vcf.gz | grep -v "^#" | \
  awk '{if(length($4)>length($5)) print}' | wc -l

echo "Variants by chromosome:"
gunzip -c clair3_output/merge_output.vcf.gz | grep -v "^#" | \
  awk '{print $1}' | sort | uniq -c
