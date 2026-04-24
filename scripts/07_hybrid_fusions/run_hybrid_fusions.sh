#!/bin/bash
set -euo pipefail

############################################################
# 07. Identification of Allelic and Non-Allelic Fusion
#     Transcripts in Hybrid Rice Using the Trans-Mate Pair Strategy
#
# Purpose:
# Identify candidate allelic and non-allelic fusion transcripts
# in the hybrid rice line SY63 using paired-end short-read RNA-seq data.
#
# Description:
# This workflow classifies read pairs according to their parental
# mapping patterns against the MH63 and ZS97 genomes, and uses
# informative cis- and trans-mate pairs to infer candidate hybrid
# fusion genes.
#
# Related figures:
# Fig. S18a-g
############################################################


############################
# A. Phase reads and assess data quality
############################

echo "Step A1: Phasing RNA-seq reads with PP2PG..."

phasing.py \
--pt=RNAseq \
--g1=toMH63.fs.bam \
--g2=toZS97.fs.bam \
--snp=MH63vsZS97.snp \
--gop1=MH63 \
--gop2=ZS97 \
--st=5

echo "Step A2: Sorting phased BAM files..."

samtools sort -@ 10 mh.bam > mh.sort.bam
samtools index mh.sort.bam

samtools sort -@ 10 zs.bam > zs.sort.bam
samtools index zs.sort.bam

echo "Step A3: Summarizing BAM correlations..."

multiBamSummary bins \
-b *sort.bam \
-o rna.multiBamSummary \
-bs 100000 \
--outRawCounts rna.multiBamSummary.rawCounts \
--scalingFactors rna.multiBamSummary.scalingFactors \
-p 30

echo "Step A4: Generating coverage tracks..."

bamCoverage \
-b mh.sort.bam \
-o mh.sort.bw \
-of bigwig \
-bs 1 \
-p 10


############################
# B. Extract read pairs with differential parental matching patterns
############################

echo "Step B1: Name-sorting parental alignment BAM files..."

samtools view -@ 4 -b -F 256 toMH63.bam | \
samtools sort -@ 4 -n - > toMH63.nameSorted.bam

samtools view -@ 4 -b -F 256 toZS97.bam | \
samtools sort -@ 4 -n - > toZS97.nameSorted.bam

echo "Step B2: Converting BAM files to BEDPE format..."

bamToBed -bedpe -ed -mate1 -i toMH63.nameSorted.bam > toMH63.nameSorted.bedpe
bamToBed -bedpe -ed -mate1 -i toZS97.nameSorted.bam > toZS97.nameSorted.bedpe

echo "Step B3: Extracting NM-tag mismatch statistics..."

samtools view toMH63.nameSorted.bam | \
perl scripts/07_hybrid_fusions/parse_bam_for_NM.pl - > toMH63.nameSorted.stats

samtools view toZS97.nameSorted.bam | \
perl scripts/07_hybrid_fusions/parse_bam_for_NM.pl - > toZS97.nameSorted.stats


############################
# C. Merge parental mismatch statistics
############################

echo "Step C1: Merging parental NM-tag summaries..."

paste \
<(cut -f 1,4,8 toMH63.nameSorted.stats) \
<(cut -f 4,8 toZS97.nameSorted.stats) \
> sy63.tagNM

echo "Step C2: Converting NM tags to parental match categories..."

cat sy63.tagNM | \
perl scripts/07_hybrid_fusions/nm_to_match.pl - > sy63.parentalAlign

echo "Step C3: Summarizing category proportions..."

cat sy63.parentalAlign | \
perl scripts/07_hybrid_fusions/make_category.pl -


############################
# D. Extract cis- and trans-mate pair reads
############################

echo "Step D1: Labeling cis/trans read-pair categories..."

cat sy63.parentalAlign | \
perl scripts/07_hybrid_fusions/get_category.pl - > sy63.cistrans

echo "Step D2: Extracting informative BEDPE read pairs..."

perl scripts/07_hybrid_fusions/get_bedpe.pl \
sy63.cistrans toMH63.nameSorted.bedpe > sy63.mh.bedpe.filter

perl scripts/07_hybrid_fusions/get_bedpe.pl \
sy63.cistrans toZS97.nameSorted.bedpe > sy63.zs.bedpe.filter

echo "Step D3: Converting BEDPE files to BED format..."

perl scripts/07_hybrid_fusions/bed12Tobed6_mh.pl \
sy63.cistrans sy63.mh.bedpe.filter > sy63.mh.bedpe.filter.bed

perl scripts/07_hybrid_fusions/bed12Tobed6_zs.pl \
sy63.cistrans sy63.zs.bedpe.filter > sy63.zs.bedpe.filter.bed


############################
# E. Identify candidate allelic and non-allelic fusion genes
############################

echo "Step E1: Retaining the longest isoform annotations..."

agat_sp_keep_longest_isoform.pl \
-f Omh63.maker.gff \
-o Omh63.longest.gff

agat_sp_keep_longest_isoform.pl \
-f Ozs97.maker.gff \
-o Ozs97.longest.gff

echo "Step E2: Adding exon order information..."

perl scripts/07_hybrid_fusions/add_exon_order.pl \
Omh63.longest.gff | \
sort -k1,1 -k4,4n -k5,5 > Omh63.exon.gff

perl scripts/07_hybrid_fusions/add_exon_order.pl \
Ozs97.longest.gff | \
sort -k1,1 -k4,4n -k5,5 > Ozs97.exon.gff

echo "Step E3: Intersecting candidate reads with exon annotations..."

intersectBed \
-a sy63.mh.bedpe.filter.bed \
-b Omh63.exon.gff \
-wao > sy63.mh.bedpe.filter.bed.exon

intersectBed \
-a sy63.zs.bedpe.filter.bed \
-b Ozs97.exon.gff \
-wao > sy63.zs.bedpe.filter.bed.exon

echo "Step E4: Determining overlapping genes..."

perl scripts/07_hybrid_fusions/reads_overlap_exon.pl \
Omh63.gene2trans.map \
sy63.mh.bedpe.filter.bed.exon > sy63.mh.readStatus

perl scripts/07_hybrid_fusions/reads_overlap_exon.pl \
Ozs97.gene2trans.map \
sy63.zs.bedpe.filter.bed.exon > sy63.zs.readStatus

echo "Step E5: Identifying candidate allelic fusion genes..."

cat sy63.mh.readStatus sy63.zs.readStatus | \
perl scripts/07_hybrid_fusions/cis_category.pl - | \
grep 'Diff_exon' > sy63.merged.cis.diffexon

cat sy63.merged.cis.diffexon | \
perl scripts/07_hybrid_fusions/cis_map_gene.pl - | \
sort -k1,1 -k2,2 > sy63.merged.cis.diffexon.exonLevel

cat sy63.merged.cis.diffexon.exonLevel | \
perl scripts/07_hybrid_fusions/convert_cis_exon_to_gene.pl - | \
sort -k1,1 > sy63.merged.cis.diffexon.geneLevel

echo "Step E6: Identifying candidate non-allelic fusion genes..."

cat sy63.mh.readStatus sy63.zs.readStatus | \
awk '$2~"trans"' | \
perl scripts/07_hybrid_fusions/merge_trans_pair.pl - > sy63.merged.trans

cat sy63.merged.trans | \
perl scripts/07_hybrid_fusions/trans_categroy.pl - | \
awk '$1=="Pass"' | \
cs-th cut -f -1 | \
perl scripts/07_hybrid_fusions/trans_find_properPair.pl - > sy63.trans.diffexon

cat sy63.trans.diffexon | \
grep -v 'multie' | \
perl scripts/07_hybrid_fusions/trans_map_gene.pl - > sy63.trans.diffexon.exonLevel

cat sy63.trans.diffexon.exonLevel | \
perl scripts/07_hybrid_fusions/convert_trans_exon_to_gene.pl - | \
sort -k1,1 > sy63.trans.diffexon.geneLevel

echo "Allelic candidates: sy63.merged.cis.diffexon.geneLevel"
echo "Non-allelic candidates: sy63.trans.diffexon.geneLevel"