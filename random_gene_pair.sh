#########################################################################
# File Name: random_gene_pair.sh
# Author: Zhu Xiaobu
# mail: z724@qq.com
# Created Time: Thu 10 Apr 2025 08:50:19 PM CST
#########################################################################
#!/bin/bash


exp_gene=$1
fusions=$2
out=$3

for n in `seq 1 5000`
do
	perl random_gene_pair.pl $1 $2 | awk -v n=$n '{x=sprintf("%04d", n); print $0"\tRand_"x}' >> $out.random_expressed_pair
done
