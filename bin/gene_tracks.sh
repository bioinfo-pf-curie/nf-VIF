#!/bin/bash

hpvlist=$1
hpvgenes=$2
prefix=$3

while IFS= read -r hpv; do
    awk -v strain=${hpv} 'NR==1{for (i=1;i<=NF;i++){header[i]=$i}}\
$1==strain{OFS="\t";for (i=1;i<=NF;i++){if($i~/\.\./){split($i,out,/\.\./); print strain,out[1],out[2],header[i]}}}' ${hpvgenes} > genes_${hpv}.bed

    awk -v strain=${hpv} '$1==strain{OFS="\t"; print $1,$2}' ${hpvgenes} > ${hpv}.size
    while IFS= read -r gene; do
	echo ${gene} | awk '{OFS="\t"; print $1,$2,$3}' > gene.coord
	genename=$(echo ${gene} | awk '{print $4}')
	genomeCoverageBed -dz -i gene.coord -g ${hpv}.size | awk '$3>0{$3=-$3*5}{OFS="\t";print $2,$3}' > ${prefix}_${genename}_${hpv}_bkp.mqc
    done < "genes_${hpv}.bed"

done < "$hpvlist"
