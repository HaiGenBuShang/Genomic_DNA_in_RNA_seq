#!/bin/bash

data_dir=${2}
file_list=${1}
ref_dir=${3}
res_dir=${4}

#Trimmomatic file
trimmo=${5}
adapter=${6}


cd $data_dir
for i in $(cat $file_list)
do

java -Xmx4g -jar ${trimmo} PE -threads 4 ${i}_1.fastq.gz ${i}_2.fastq.gz ${i}_1.trimmed.fastq.gz ${i}_1.unpaired.fastq.gz ${i}_2.trimmed.fastq.gz ${i}_2.unpaired.fastq.gz ILLUMINACLIP:${adapter}:2:30:10 LEADING:10 TRAILING:10 HEADCROP:10 SLIDINGWINDOW:4:15 MINLEN:36 

mkdir -p ${7}/${i}
mkdir -p ${res_dir}/${i}
sam_res_dir=${res_dir}/${i}
bam_dir=${7}/${i}
hisat2 -q -p 4 -x ${ref_dir}/genome -1 ${i}_1.trimmed.fastq.gz -2 ${i}_2.trimmed.fastq.gz -S ${bam_dir}/${i}.sam

samtools view -b ${bam_dir}/${i}.sam >  ${bam_dir}/${i}.bam
rm ${bam_dir}/${i}.sam
samtools sort -@ 3 ${bam_dir}/${i}.bam > ${bam_dir}/${i}.sorted.bam
#samtools rmdup -S ${bam_dir}/${i}.sorted.bam ${bam_dir}/${i}.sorted.rmdup.bam
rm ${bam_dir}/${i}.bam
done
echo "gDNA_mapping.sh finished! "
