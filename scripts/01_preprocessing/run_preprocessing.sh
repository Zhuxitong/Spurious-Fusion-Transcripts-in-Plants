#!/bin/bash
set -euo pipefail

############################################################
# 01. Sequencing Data Preprocessing
# Purpose:
# Preprocess PacBio Iso-Seq and Nanopore long-read datasets
# for downstream fusion transcript analyses.
#
# Related figure:
# Supplementary Figure S1
############################################################

############################
# A. PacBio Iso-Seq
############################

# Input:
#   isoseq.subreads.bam
#   isoseq.subreads.bam.pbi
#   isoseq.subreads.bam.xml
#   test_data/primers.fasta

echo "Step A1: Generate CCS reads..."
ccs isoseq.subreads.bam isoseq.ccs.bam \
-j 3 \
--min-length 50 \
--min-rq 0.9 \
--report-file isoseq.ccs.report

echo "Step A2: Primer removal and full-length read detection..."
lima isoseq.ccs.bam test_data/primers.fasta isoseq.fl.bam \
--isoseq -j 5

echo "Step A3: Generate FLNC reads..."
isoseq refine \
isoseq.fl.primer_5p--primer_3p.bam \
test_data/primers.fasta \
isoseq.flnc.bam \
--require-polya -j 5

echo "Step A4: Convert FLNC BAM to FASTA..."
samtools fasta isoseq.flnc.bam > isoseq.flnc.fasta


############################
# B. Nanopore direct RNA
############################

# Input:
#   dRNA.fastq.gz
#   short_reads.fq.gz
#   Omh63.genomic.fa

echo "Step B1: Generate read statistics..."
seqkit stats -b dRNA.fastq.gz

echo "Step B2: Quality filtering..."
chopper -l 150 -q 7 -i dRNA.fastq.gz | pigz - > dRNA.qc.fastq.gz

echo "Step B3: Build short-read index for FMLRC2..."
pigz -c -d short_reads.fq.gz | \
awk 'NR % 4 == 2' | \
sort | \
tr NT TN | \
ropebwt2 -LR | \
tr NT TN | \
msbwt2-convert short_reads_fmlrc2.npy

echo "Step B4: Correct Nanopore reads..."
fmlrc2 -C 10 -t 2 \
short_reads_fmlrc2.npy \
dRNA.qc.fastq.gz \
dRNA.correct.fasta

echo "Step B5: Align corrected reads to genome..."
minimap2 -ax splice -uf -k14 -t 2 Omh63.genomic.fa dRNA.correct.fasta | \
samtools sort --threads 2 - | \
samtools view --threads 2 -b - > dRNA.2mh.bam