#!/bin/bash
set -euo pipefail

############################################################
# 08. Additional Analyses
#
# This script includes:
#   A. Overlap between fusion transcripts and 3D interaction data
#   B. Assembly of the MH63 primary genome using HiFi reads
#   C. Comparison of different genome annotation versions
############################################################


############################
# A. Overlap with 3D interaction data
############################

echo "Step A1: Intersecting fusion transcript pairs with 3D interaction data..."

pairToPair \
-a test_data/isoseq.fusions.promoter.bedpe \
-b test_data/dna2rna_H3K4me3_rs3.bedpe \
-is > isoseq.dna2rna_H3K4me3.observed

echo "Step A2: Generating matched random gene pairs and calculating random overlaps..."

sh random_gene_pair.sh

pairToPair \
-a test_data/random_expressed_pair.promoter.bedpe \
-b test_data/dna2rna_H3K4me3_rs3.bedpe \
-is > dna2rna_H3K4me3.random


############################
# B. MH63 HiFi assembly
############################

echo "Step B1: Assembling the MH63 primary genome using HiFi reads..."

hifiasm -o MH63.asm -t 32 SRR10238608.fastq.gz


############################
# C. Comparison of annotation versions
############################

echo "Step C1: Calculating overall annotation statistics..."

agat_sp_statistics.pl \
--gff Omh63.gff3 \
-g Omh63.genomic.fa \
-o mh63.stats

agat_sp_statistics.pl \
--gff Omh63_man.gff3 \
-g Omh63.genomic.fa \
-o mh63_man.stats


echo "Step C2: Converting exon annotations to BED format..."

agat_convert_sp_gff2bed.pl \
--gff Omh63.gff3 \
--sub exon \
-o Omh63.bed

agat_convert_sp_gff2bed.pl \
--gff Omh63_man.gff3 \
--sub exon \
-o Omh63_man.bed


echo "Step C3: Running ParsEval to compare annotation structures..."

parseval -d -o parseval.out Omh63.bed Omh63_man.bed &> parseval.log


echo "Step C4: Identifying split-type gene annotations..."

intersectBed \
-a <(awk '$3=="gene"' Omh63.gff3) \
-b <(awk '$3=="gene"' Omh63_man.gff3) \
-wao -s | \
sed 's/ID=//g' | \
awk '$18!="."' | \
csvtk_tH fold -f 9 -v 18 -s ";" | \
awk '{n=sub(/;/,";",$2); split($2,a,";"); if(n>=1){print $0"\t"length(a)}}' \
> v1_split_in_v2.id


echo "Step C5: Identifying fused-type gene annotations..."

intersectBed \
-a <(awk '$3=="gene"' Omh63_man.gff3) \
-b <(awk '$3=="gene"' Omh63.gff3) \
-wao -s | \
sed 's/ID=//g' | \
awk '$18!="."' | \
csvtk_tH fold -f 9 -v 18 -s ";" | \
awk '{n=sub(/;/,";",$2); split($2,a,";"); if(n>=1){print $0"\t"length(a)}}' \
> v1_fuse_in_v2.id


echo "Step C6: Summarizing new transcripts and detailed locus comparisons..."

perl resolve_parseval_out.newTranscripts.pl parseval.out | \
awk '$NF!=0' > parseval.out.newTrans

perl resolve_parseval_out.pl parseval.out > parseval.out.locusNum

perl resolve_parseval_out.comparison.pl parseval.out > parseval.out.comparison