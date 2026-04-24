#!/bin/bash
set -euo pipefail

############################################################
# 04. Classification of Fusion Breakpoint Types
#
# Purpose:
# Classify putative fusion breakpoints into four categories:
# SHS, Joined, Unknown, and Misaligned.
#
# Related figures:
# Fig. 2a-b, Fig. S4c, Fig. S5a
############################################################


############################
# A. Convert breakpoint coordinates
############################

echo "Step A1: Extract genomic breakpoint coordinates..."

cut -f 1,2,3,4,6,7 mh63_flagleaf.pc | \
sed "s/:/\t/g" | \
perl scripts/04_breakpoint_classification/get_breakpoints_genomic_coords.pl - | \
sortBed -i - > mh63_flagleaf.pc.breakpoints.genomic.bed

echo "Step A2: Please run the accompanying R script (convert.R) to map genomic breakpoints to transcript coordinates."
echo "Expected output: mh63_flagleaf.pc.breakpoints.trans.txt"

echo "Step A3: Convert transcript-space breakpoints to BED format..."

sed '1d' mh63_flagleaf.pc.breakpoints.trans.txt | \
awk '{print $2"\t"$3-1"\t"$4"\t"$9}' | \
sortBed > mh63_flagleaf.pc.breakpoints.trans.bed


############################
# B. Prepare breakpoint and coverage files
############################

echo "Step B1: Extract breakpoint coordinates on fusion reads..."

perl scripts/04_breakpoint_classification/get_breakpoints_on_read.pl \
mh63_flagleaf.pc Isoseq.csv > mh63_flagleaf.pc.reads.breakpoints.bed

echo "Step B2: Merge parental-gene alignments..."

mergeBed -d -30 -c 4 -o distinct \
-i mh63_flagleaf.pc.reads.gene.bed > mh63_flagleaf.pc.reads.gene.merged.bed

echo "Step B3: Build FASTA index and length table..."

samtools faidx mh63_flagleaf.pc.reads.fa
cut -f 1,2 mh63_flagleaf.pc.reads.fa.fai > mh63_flagleaf.pc.reads.fa.length


############################
# C. Compute coverage and classify breakpoints
############################

echo "Step C1: Compute parental-gene coverage..."

genomeCoverageBed \
-i mh63_flagleaf.pc.reads.gene.merged.bed \
-g mh63_flagleaf.pc.reads.fa.length \
-bga > mh63_flagleaf.pc.reads.gene.coverage

echo "Step C2: Identify SHS breakpoints with overlap >3 bp..."

intersectBed \
-a <(cut -f 1-3 mh63_flagleaf.pc.reads.breakpoints.bed) \
-b mh63_flagleaf.pc.reads.gene.coverage \
-wao | \
awk '$7==2{print $1"\t"$5"\t"$6"\t"$6-$5}' | \
awk '$4>3' > mh63_flagleaf.pc.reads.out.shsGreat3.bed

cut -f 1 mh63_flagleaf.pc.reads.out.shsGreat3.bed > mh63_flagleaf.pc.reads.out.shsGreat3.id

intersectBed \
-a <(cut -f 1-3 mh63_flagleaf.pc.reads.breakpoints.bed) \
-b mh63_flagleaf.pc.reads.gene.coverage \
-wao | \
awk '$7==2{print $1"\t"$5"\t"$6"\t"$6-$5}' \
> mh63_flagleaf.pc.reads.out.interBed

awk '$4<=3' mh63_flagleaf.pc.reads.out.interBed > mh63_flagleaf.pc.reads.out.join1.bed
cut -f 1 mh63_flagleaf.pc.reads.out.interBed > tmp.id

awk '$4==0' mh63_flagleaf.pc.reads.gene.coverage > mh63_flagleaf.pc.reads.out.coverage0
awk '$4==2' mh63_flagleaf.pc.reads.gene.coverage > mh63_flagleaf.pc.reads.out.coverage2


echo "Step C3: Identify Joined breakpoints..."

intersectBed \
-a <(cut -f 1-3 mh63_flagleaf.pc.reads.breakpoints.bed) \
-b mh63_flagleaf.pc.reads.gene.coverage \
-wao | \
awk '$7<2' | \
cut -f 1-3 | \
cs-th grep -f 1 -v -P tmp.id | \
sort | uniq | sortBed | \
closestBed -a - -b mh63_flagleaf.pc.reads.out.coverage0 -d | \
awk '$4=="."' | \
cut -f 1-3 > mh63_flagleaf.pc.reads.out.join2.bed

intersectBed \
-a <(cut -f 1-3 mh63_flagleaf.pc.reads.breakpoints.bed) \
-b mh63_flagleaf.pc.reads.gene.coverage \
-wao | \
awk '$7<2' | \
cut -f 1-3 | \
cs-th grep -f 1 -v -P tmp.id | \
sort | uniq | sortBed | \
closestBed -a - -b mh63_flagleaf.pc.reads.out.coverage0 -d | \
awk '$4!="."' | \
awk '{print $1"\t"$2"\t"$3"\t"$6-$5"\t"$8}' | \
awk '$5>=30' | \
closestBed -a - -b mh63_flagleaf.pc.reads.out.coverage2 -d | \
awk '{print $1"\t"$2"\t"$3"\t"$8-$7"\t"$10}' | \
awk '$4<=3 || $5>10' > mh63_flagleaf.pc.reads.out.join3.bed

intersectBed \
-a <(cut -f 1-3 mh63_flagleaf.pc.reads.breakpoints.bed) \
-b mh63_flagleaf.pc.reads.gene.coverage \
-wao | \
cut -f 1-3 | \
cs-th grep -f 1 -v -P tmp.id | \
sort | uniq | sortBed | \
closestBed -a - -b mh63_flagleaf.pc.reads.out.coverage0 -d | \
awk '$4!="."' | \
awk '{print $1"\t"$2"\t"$3"\t"$6-$5"\t"$8}' | \
awk '$5<30 && $4<30' | \
awk '$4==1 || $4==2' > mh63_flagleaf.pc.reads.out.join4.bed


echo "Step C4: Identify additional SHS, Unknown, and Misaligned breakpoints..."

intersectBed \
-a <(cut -f 1-3 mh63_flagleaf.pc.reads.breakpoints.bed) \
-b mh63_flagleaf.pc.reads.gene.coverage \
-wao | \
cut -f 1-3 | \
cs-th grep -f 1 -v -P tmp.id | \
sort | uniq | sortBed | \
closestBed -a - -b mh63_flagleaf.pc.reads.out.coverage0 -d | \
awk '$4!="."' | \
awk '{print $1"\t"$2"\t"$3"\t"$6-$5"\t"$8}' | \
awk '$5>=30' | \
closestBed -a - -b mh63_flagleaf.pc.reads.out.coverage2 -d | \
awk '{print $1"\t"$7"\t"$8"\t"$8-$7"\t"$10}' | \
awk '$4>3 && $5<=10' > mh63_flagleaf.pc.reads.out.shsGreat3_1.bed

intersectBed \
-a <(cut -f 1-3 mh63_flagleaf.pc.reads.breakpoints.bed) \
-b mh63_flagleaf.pc.reads.gene.coverage \
-wao | \
cut -f 1-3 | \
cs-th grep -f 1 -v -P tmp.id | \
sort | uniq | sortBed | \
closestBed -a - -b mh63_flagleaf.pc.reads.out.coverage0 -d | \
awk '$4!="."' | \
awk '{print $1"\t"$2"\t"$3"\t"$6-$5"\t"$8}' | \
awk '$5<30 && $4<30' | \
awk '$4!=1 && $4!=2' > mh63_flagleaf.pc.reads.out.unknown.bed

intersectBed \
-a <(cut -f 1-3 mh63_flagleaf.pc.reads.breakpoints.bed) \
-b mh63_flagleaf.pc.reads.gene.coverage \
-wao | \
cut -f 1-3 | \
cs-th grep -f 1 -v -P tmp.id | \
sort | uniq | sortBed | \
closestBed -a - -b mh63_flagleaf.pc.reads.out.coverage0 -d | \
awk '$4!="."' | \
awk '{print $1"\t"$2"\t"$3"\t"$6-$5"\t"$8}' | \
awk '$4>=30 && $5<30' > mh63_flagleaf.pc.reads.out.error.bed


############################
# D. Merge final outputs
############################

echo "Step D1: Merge SHS breakpoints..."

cat mh63_flagleaf.pc.reads.out.shsGreat3.bed \
mh63_flagleaf.pc.reads.out.shsGreat3_1.bed | \
cut -f 1-4 > mh63_flagleaf.pc.reads.out.shs.bed

echo "Step D2: Merge Joined breakpoints..."

cat mh63_flagleaf.pc.reads.out.join?.bed | \
cut -f 1-3 | \
awk '{print $0"\t0"}' > mh63_flagleaf.pc.reads.out.join.bed