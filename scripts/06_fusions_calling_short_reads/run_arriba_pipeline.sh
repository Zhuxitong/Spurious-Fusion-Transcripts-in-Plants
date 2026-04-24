#!/bin/bash
set -euo pipefail

############################################################
# 06. Fusion Transcript Calling from Short-Read RNA-seq Data
#
# Purpose:
# Detect putative fusion transcripts using STAR + Arriba.
# Used for adjacent-gene fusion analyses and comparison
# with long-read sequencing results.
#
# Related figures:
# Fig. 5a, Fig. S11, Fig. S13a
#
# Data accession:
# PRJNA1291274
############################################################


############################
# A. Build STAR genome index
############################

echo "Step A1: Building STAR genome index..."

STAR \
--runMode genomeGenerate \
--genomeDir ./ \
--genomeFastaFiles Omh63.fa \
--sjdbGTFfile Omh63.gtf \
--runThreadN 4 \
--genomeSAindexNbases 13


############################
# B. Align short reads
############################

echo "Step B1: Aligning paired-end RNA-seq reads with STAR..."

STAR \
--runMode alignReads \
--runThreadN 5 \
--genomeDir mh63 \
--genomeLoad NoSharedMemory \
--readFilesIn data1.fq.gz data2.fq.gz \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outFilterMultimapNmax 50 \
--peOverlapNbasesMin 10 \
--alignSplicedMateMapLminOverLmate 0.5 \
--alignSJstitchMismatchNmax 5 -1 5 5 \
--chimSegmentMin 10 \
--chimOutType WithinBAM HardClip \
--chimJunctionOverhangMin 10 \
--chimScoreDropMax 30 \
--chimScoreJunctionNonGTAG 0 \
--chimScoreSeparation 1 \
--chimSegmentReadGapMax 3 \
--chimMultimapNmax 50 \
--outFileNamePrefix output_pre


############################
# C. Run Arriba
############################

echo "Step C1: Detecting fusion transcripts using Arriba..."

arriba \
-x output_preAligned.sortedByCoord.out.bam \
-g Omh63.gtf \
-a Omh63.fa \
-o output_pre.arribaFusion.tsv \
-O output_pre.arribaFusion.discarded \
-f blacklist \
-i Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12


############################
# D. Merge and filter outputs
############################

echo "Step D1: Reformatting Arriba output files..."

for i in *.tsv; do
pre=$(echo $i | cut -f 1 -d '.')
perl reformat_tsv.pl $i | \
awk -v n=$pre '{print $0"\t"n}' > $pre.basic
done


echo "Step D2: Merging replicates..."

cat *.basic | perl merge_rep.pl - > flagleaf.arriba


echo "Step D3: Filtering repetitive candidates..."

for i in *.arriba; do
pre=$(basename $i | cut -f 1 -d '.')
cat $i | perl filter_repeats.pl - > $pre.arriba
done


echo "Step D4: Merging all treatment groups..."

cat flagleaf_HTLD.arriba \
flagleaf_HTSD.arriba \
flagleaf_LTLD.arriba \
flagleaf_LTSD.arriba | \
perl merge_flagleaf.pl - | \
csvtk -tH sort -k 13:u -L 13:ranks > flagleaf.result.filter1_rep