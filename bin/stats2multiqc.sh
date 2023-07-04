#!/bin/bash

splan=$1

## Catch sample names
all_samples=$(awk -F, '{print $1}' $splan)

echo -e "sample_id,sample_name,reads_number,reads_with_adapter,percent_with_adapter,reads_on_ctrl,percent_on_ctrl,unique_HPV,percent_unique_HPV,multi_HPV,percent_multi_HPV,genotypes,status" > mq_main.stats

for sample in $all_samples
do

    ##id
    sname=$(awk -F, -v sname=$sample '$1==sname{print $2}' $splan)

    ##n_reads
    trimFile=$(find trimming/ -name "${sample}*" -name '*trimming_report.txt' | awk 'NR==1{print}')
    n_reads=$(grep "Total reads processed:" ${trimFile} | cut -d: -f 2 | cut -d"(" -f 1 | sed -e 's/,//g' -e 's/ //g')
    
    ##n_reads_trimmed
    n_with_adapt=$(grep "Reads with adapters:" ${trimFile} | cut -d: -f 2 | cut -d"(" -f 1 | sed -e 's/,//g' -e 's/ //g')
    p_with_adapt=$(echo "${n_with_adapt} ${n_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')

    ##n_reads_ctrl
    n_ctrl=$(awk '$1!="unmapped"{s=s+$2}END{print s}' ctrl/${sample}_ctrl.stats)
    p_ctrl=$(echo "${n_ctrl} ${n_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')

    ##n_reads_hpvs
    n_unique_hpv=$(grep "aligned concordantly exactly 1 time" hpv/${sample}_hpvs_bowtie2.log | cut -d"(" -f 1 | sed -e 's/ //g')
    n_multi_hpv=$(grep "aligned concordantly >1 time" hpv/${sample}_hpvs_bowtie2.log | cut -d"(" -f 1 | sed -e 's/ //g')
    p_unique_hpv=$(echo "${n_unique_hpv} ${n_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    p_multi_hpv=$(echo "${n_multi_hpv} ${n_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    p_hpv=$(echo "$p_unique_hpv + $p_multi_hpv" | bc)

    ##detected_genotypes
    geno=$(cut -d, -f2 hpv/${sample}_HPVgenotyping.filtered | tr '\n' '|' | sed -e 's/|$//')
    
    ##n_breakpoints
    n_hq_bkp=-1
    for bkpFile in $(find hpv/ -name "${sample}*" -name '*bkptable_filtered.csv')
    do
	if [[ $n_hq_bkp == -1 ]];then
	    n_hq_bkp=0
	fi
	genotype=$(basename $bkpFile | sed -e 's/'${sample}'-//' -e 's/_bkptable_filtered.csv//')
	#nb_geno=$(awk -F"\t" -v geno=${genotype} '$1==geno{print $2}' hpv/${sample}_HPVgenotyping.stats)
	nb_geno=$(awk -F',' -v geno=${genotype} '$3==geno{print $4}' hpv/${sample}-${genotype}_coverage.stats)
	min_count=$(echo $nb_geno | awk -v minFrac=0.0001 '{printf "%d", $nb_gen*minFrac}')
	nbbkp=$(awk -F, -v minCount=${min_count} '$8==10 && $7>minCount{print}' $bkpFile | wc -l)
	n_hq_bkp=$(($n_hq_bkp + $nbbkp))
	echo "$genotype - $nb_geno - $min_count - $nbbkp"
    done
    
    ##status
    if [[ ${n_hq_bkp} == -1 ]]; then
	status="UNKNOWN"
    elif [[ $(echo "${p_unique_hpv} >= 1" |bc -l) == 1 && $(echo "${n_hq_bkp} >= 1" |bc -l) == 1 ]];then
	status="INTEGRATED"
    elif [[ $(echo "${p_unique_hpv} >= 1" |bc -l) == 1 ]]; then
	status="EPI"
    else
	status="NEG"
    fi

    echo -e ${sample},${sname},${n_reads},${n_with_adapt},${p_with_adapt},${n_ctrl},${p_ctrl},${n_unique_hpv},${p_unique_hpv},${n_multi_hpv},${p_multi_hpv},${geno},${status} >> mq_main.stats
done
