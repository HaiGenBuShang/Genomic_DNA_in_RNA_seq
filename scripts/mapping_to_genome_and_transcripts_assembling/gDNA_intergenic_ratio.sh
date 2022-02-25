#!/bin/bash

data_dir=${2}
file_list=${1}
ref_dir=${3}
res_dir=${4}

#Trimmomatic file
trimmo=${5}
adapter=${6}

#annotation file
#featureCounts count GeneID column as counting meta-feature
anno=${7}
intergenic=${8}

cd $data_dir
for i in $(cat $file_list)
do

sam_res_dir=${res_dir}/${i}
bam_dir=${9}/${i}

#to report the mapped reads to the whole genome and uniquely mapped intergenic region
#samtools stats ${bam_dir}/${i}.sorted.rmdup.bam > ${sam_res_dir}/${i}.sorted.rmdup.bam.stats.txt
#samtools stats ${sam_res_dir}/${i}.sorted.rmdup.intergenic.unique.bam > ${sam_res_dir}/${i}.sorted.rmdup.intergenic.unique.bam.stats.txt


#to report the mapped reads to the whole genome and uniquely mapped intergenic region
samtools stats ${bam_dir}/${i}.sorted.bam > ${sam_res_dir}/${i}.sorted.bam.stats.txt
samtools stats ${sam_res_dir}/${i}.sorted.intergenic.unique.bam > ${sam_res_dir}/${i}.sorted.intergenic.unique.bam.stats.txt

done

echo "gDNA_intergenic_ratio.sh finished! "
