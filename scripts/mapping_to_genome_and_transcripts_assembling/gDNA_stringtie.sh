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

current_dir=$(pwd)

cd $data_dir
for i in $(cat $file_list)
do

sam_res_dir=${res_dir}/${i}
bam_dir=${9}/${i}

cd ${9}
#The stringtie input should be coordinate sorted
#stringtie ${bam_dir}/${i}.sorted.rmdup.bam -p 3 -a 10 -m 200 -g 50 -c 2.5 -j 1 -l STRG -M 0.95 -G ${anno} -e -o ${bam_dir}/${i}.transcripts.stringtie.gtf -A ${bam_dir}/${i}.gene.abundance.txt -b ./BallgownTable/${i}
stringtie ${bam_dir}/${i}.sorted.bam -p 3 -a 10 -m 200 -g 50 -c 2.5 -j 1 -l STRG -M 0.95 -G ${anno} -e -o ${bam_dir}/${i}.transcripts.stringtie.gtf -A ${bam_dir}/${i}.gene.abundance.txt -b ./BallgownTable/${i}
done
echo "gDNA_stringtie.sh finished!"

printf "${current_dir}/Stringtie finished! \n" >> ${current_dir}/wait_to_delete_automatically.txt
