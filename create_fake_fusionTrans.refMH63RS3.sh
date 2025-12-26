#########################################################################
# File Name: fuckFusion.sh
# Author: Zhu Xiaobu
# mail: z724@qq.com
# Created Time: Thu 03 Nov 2022 11:16:09 PM CST
# Modified Time: 2025-02-26
#########################################################################
#!/bin/bash


# MH63 genomic info 
locusFa=/home/xtzhu/works/2023-09-15_fusionTranscripts/results/05init_jaffal_results/visual_files_mh63/Omh63_RS3.fasta
locusAnno=/home/xtzhu/works/2023-09-15_fusionTranscripts/results/05init_jaffal_results/visual_files_mh63/Omh63.transcriptAnno.gtf
locusLen=/home/xtzhu/works/2023-09-15_fusionTranscripts/results/05init_jaffal_results/visual_files_mh63/Omh63.geneLocus.fa.length
locusSNPs=/home/xtzhu/works/2023-09-15_fusionTranscripts/results/05init_jaffal_results/visual_files_mh63/Omh63.transSNPs.bed
locusGenoPos=/home/xtzhu/works/2023-09-15_fusionTranscripts/results/05init_jaffal_results/visual_files_mh63/Omh63_RS3.txt


# breakpoint on trans
## need to generate from jaffal results. Which means must be changed when the input fusion events are changed.
fusionMapList=/home/xtzhu/works/2023-09-15_fusionTranscripts/results/05init_jaffal_results/02drs_jaffal/mh63_drs.pc.breakpoints.trans.bed


#fusionEvent='OsMH_01G0242900:OsMH_06G0262900#Chr01#Chr06'
fusionEvent=$1
fusionGene=$(echo ${fusionEvent} | cut -f 1 -d '#' | sed 's/:/#/')
parentGene1=$(echo ${fusionGene} | cut -f 1 -d '#')
parentGene2=$(echo ${fusionGene} | cut -f 2 -d '#')
chr_1=$(echo ${fusionEvent} | cut -f 2 -d '#')
chr_2=$(echo ${fusionEvent} | cut -f 3 -d '#')
#echo -e "${parentGene1}\n${parentGene2}\n"
out=$2

# fasta
seqkit grep -p ${parentGene1} ${locusFa} > tmp_parentGene1.fa
seqkit grep -p ${parentGene2} ${locusFa} | sed "s/${parentGene2}/${parentGene1}/" > tmp_parentGene2.fa
seqkit concat tmp_parentGene1.fa tmp_parentGene2.fa | seqkit replace -p '(.+)' -r "$fusionGene" >> $out.visual.fa
rm tmp_parentGene1.fa tmp_parentGene2.fa


# anno
grep -P "${parentGene1}\t" ${locusAnno} | awk -v id=${fusionGene} -v FS="\t" -v OFS="\t" '{$1=id; print $0}' >> $out.visual.gtf
gene1End=$(grep -P "${parentGene1}\t" ${locusLen} | cut -f 2)
echo -e "${fusionGene}\t0\t${gene1End}\t${chr_1}\t0\t.\t0\t${gene1End}\t255,0,0"  >> $out.visual.gene_boundary.bed 
gene2end=$(grep -P "${parentGene2}\t" ${locusLen} | cut -f 2)
awk -v id=${fusionGene} -v pos1=${gene1End} -v pos2=${gene2end} -v chr=${chr_2} 'BEGIN{print id"\t"pos1"\t"pos1+pos2"\t"chr"\t1\t.\t"pos1"\t"pos1+pos2"\t0,186,255"}' >> $out.visual.gene_boundary.bed
grep -P "${parentGene2}\t" ${locusAnno} | awk -v id=${fusionGene} -v pos=${gene1End} -v FS="\t" -v OFS="\t" '{$1=id; $4=$4+pos; $5=$5+pos; print $0}' >> $out.visual.gtf


# SNPs
grep -P "${parentGene1}\t" ${locusSNPs} | awk -v id=${fusionGene} -v FS="\t" -v OFS="\t" '{$1=id; print $0"\t0\t.\t"$2"\t"$3"\t120,40,249"}' >> $out.visual.SNPs.bed
grep -P "${parentGene2}\t" ${locusSNPs} | awk -v id=${fusionGene} -v pos=${gene1End} -v FS="\t" -v OFS="\t" '{$1=id; $2=$2+pos; $3=$3+pos; print $0"\t1\t.\t"$2"\t"$3"\t105,179,243"}' >> $out.visual.SNPs.bed


# predicted breakPoint
grep -P "${fusionGene}" ${fusionMapList} | grep -P "${parentGene1}\t" | awk -v id=${fusionGene} -v FS="\t" -v OFS="\t" '{$1=id; print $0}' >> $out.visual.breakpoint.bed
grep -P "${fusionGene}" ${fusionMapList} | grep -P "${parentGene2}\t" | awk -v id=${fusionGene} -v pos=${gene1End} -v FS="\t" -v OFS="\t" '{$1=id; $2=$2+pos; $3=$3+pos; print $0}' >> $out.visual.breakpoint.bed


# Any other tracks can be generated!
## Such as protein domain. Need the location on transcript.
