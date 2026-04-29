#!/bin/bash
set -euo pipefail

############################################################
# 02. Fusion Transcript Identification Using Optimized JAFFAL
#
# Purpose:
# Detect putative fusion transcripts from plant long-read data
# using a modified JAFFAL pipeline.
#
# Related figures:
# Fig. 1b-d, S2, S3, S4a-b, S14b
############################################################


############################
# A. Install JAFFAL
############################

echo "Installing optimized JAFFAL..."

unzip JAFFA-version-2.3.zip
cd JAFFA-version-2.3
./install_linux64.sh
cd ..


############################
# B. Build JAFFAL Plant Database
############################

echo "Building custom JAFFAL reference database..."

mkdir -p JAFFAL_DB_MH63RS3
cd JAFFAL_DB_MH63RS3

wget http://rice.hzau.edu.cn/rice_rs3/download_ext/MH63RS3.fasta.gz
wget http://rice.hzau.edu.cn/rice_rs3/download_ext/MH63RS3_nonTE.gff3.tar.gz

mv MH63RS3.fasta Omh63.fa
gffread MH63RS3_nonTE.gff3 -T -o Omh63.gtf

# Generate transcript FASTA
gffread Omh63.gtf -g Omh63.fa -w Omh63_RS3.fasta

# Reformat transcript FASTA
cat Omh63_RS3.fasta | \
awk '{if(/^>/){print $1}else{print}}' | \
$JAFFAL_HOME/tools/bin/reformat \
fastawrap=0 in=stdin.fa out=stdout.fa > Omh63_RS3.fa

# GenePred table
gtfToGenePred Omh63.gtf Omh63.gpd -genePredExt

awk '{print "1\t"$0}' Omh63.gpd | \
csvtk -tH add-header \
-n '#bin,name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds,score,name2,cdsStartStat,cdsEndStat,exonFrames' \
> Omh63_RS3.tab

# Exon BED
grep -P "\texon\t" Omh63.gtf | \
awk '{print $1"\t"$4-1"\t"$5"\t"$10"\t.\t"$7}' | \
sortBed > Omh63_RS3.bed

# Mask genome
bedtools maskfasta \
-fi Omh63.fa \
-fo Masked_Omh63.fa \
-bed Omh63_RS3.bed

# Build indices
bowtie2-build Omh63_RS3.fa Omh63_RS3
bowtie2-build Masked_Omh63.fa Masked_Omh63

makeblastdb \
-in Omh63_RS3.fa \
-dbtype nucl \
-out Omh63_RS3_blast

touch known_fusions.txt

cd ..


############################
# C. Validate JAFFAL Database
############################

echo "Testing database with simulated fusion transcripts..."

~/PATH_to_FusionTranscriptSimulator/FusionSimulatorToolkit-master/FusionTranscriptSimulator \
MH63RS3_nonTE.gff3 Omh63.fa 20 > test.fasta

$JAFFAL_HOME/tools/bin/bpipe run -n 1 \
-p genome=Omh63 \
-p annotation=RS3 \
-p outputName=test \
-p refBase=~/PATH_to_JAFFAL_DB/JAFFAL_DB_MH63RS3 \
$JAFFAL_HOME/JAFFA_direct.groovy \
test.fasta


############################
# D. Run JAFFAL on Real Data
############################

echo "Running JAFFAL on Iso-Seq data..."

$JAFFAL_HOME/tools/bin/bpipe run -n 10 \
-p refBase=~/PATH_to_JAFFAL_DB/JAFFAL_DB_MH63RS3 \
-p genome=Omh63 \
-p annotation=RS3 \
-p outputName=Isoseq \
$JAFFAL_HOME/JAFFA_direct.groovy \
Isoseq.fasta


############################
# E. Reformat Output
############################

echo "Reformatting JAFFAL output..."

sed '1d' Isoseq.csv | \
perl scripts/02_fusion_calling/reformat_jaffal.pl - | \
csvtk -tH sort -k 5:u -L 5:ranks > mh63_flagleaf.pc

sed '1d' Isoseq.csv | \
perl scripts/02_fusion_calling/convert_jaffal_to_bedpe.pl - | \
sortBed -i - > mh63_flagleaf.pc.bedpe

awk '{print $0"\t"$10"|"$1":"$2"|"$11}' \
mh63_flagleaf.pc > mh63_flagleaf.pc.uniqueID
