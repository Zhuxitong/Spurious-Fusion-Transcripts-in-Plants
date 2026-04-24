#!/bin/bash
set -euo pipefail

############################################################
# 03. Visualization of Fusion Transcripts and Breakpoints
#
# Purpose:
# Generate IGV-compatible files for visual inspection of
# putative fusion transcripts and breakpoint structures.
#
# Related figures:
# Fig. 2c
# Fig. S6a-b, S7, S9, S12, S13b, S15, S16,
# Fig. S19b-c, S21b, S22
############################################################


############################
# A. Build pseudo-fusion IGV references
############################

echo "Preparing visualization input..."

cut -f 1-4 mh63_flagleaf.pc | \
sed "s/:/\t/g" | \
cut -f 1-3,5 | \
sed "s/\t/:/" | \
sed "s/\t/#/g" | \
sort | uniq > mh63_flagleaf.pc.visual.input


echo "Generating pseudo-fusion references..."

cat mh63_flagleaf.pc.visual.input | \
while read line; do
sh scripts/03_fusion_and_breakpoints_visualization/create_fake_fusionTrans.refMH63RS3.sh \
$line mh63_flagleaf.pc
done


############################
# Align fusion-supporting reads
############################

echo "Aligning JAFFAL fusion reads..."

cp ~/PATH_to_JAFFAL_results/Isoseq.fusions.fa \
mh63_flagleaf.jaffal.fusions.fa

minimap2 -ax splice:hq -uf \
mh63_flagleaf.pc.visual.fa \
mh63_flagleaf.jaffal.fusions.fa | \
samtools view -b | \
samtools sort > mh63_flagleaf.jaffal.fusions.bam

samtools index mh63_flagleaf.jaffal.fusions.bam


############################
# B. Breakpoint structure visualization
############################

echo "Preparing fusion transcript sequences..."

awk '{print $10"|"$1":"$2"|"$11}' \
mh63_flagleaf.pc > mh63_flagleaf.pc.reads.id

cat Isoseq.fasta | \
fasta_formatter -w 80 > mh63_flagleaf.pc.reads.fa


echo "Extracting parental gene sequences..."

cut -f 1,2 mh63_flagleaf.pc | \
sed "s/\t/\n/" | \
sort | uniq > mh63_flagleaf.pc.reads.geneID

seqtk subseq Omh63.geneLocus.fa \
mh63_flagleaf.pc.reads.geneID | \
fasta_formatter -w 80 > mh63_flagleaf.pc.reads.gene.fa


echo "Running BLAST alignment..."

blastn \
-subject mh63_flagleaf.pc.reads.fa \
-query mh63_flagleaf.pc.reads.gene.fa \
-outfmt '6 qseqid qlen qstart qend sstrand sseqid slen sstart send nident length bitscore' \
-perc_identity 95 \
> mh63_flagleaf.pc.reads.gene.blastn


echo "Generating BED tracks..."

awk '{
if($5=="minus"){s=$9-1;e=$8}
else{s=$8-1;e=$9};
print $6"\t"s"\t"e"\t"$1
}' mh63_flagleaf.pc.reads.gene.blastn | \
sortBed > mh63_flagleaf.pc.reads.gene.bed