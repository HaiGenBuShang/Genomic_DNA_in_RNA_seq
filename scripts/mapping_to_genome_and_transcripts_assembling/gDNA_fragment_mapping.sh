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
#to obtain the unique reads mapped to intergenic region
#samtools view -L ${intergenic} -h ${bam_dir}/${i}.sorted.rmdup.bam| grep -P '^\@|NH:i:1\b'| samtools view -b > ${sam_res_dir}/${i}.sorted.rmdup.intergenic.unique.bam

#read paired were mapped
#samtools view -b -F12 ${sam_res_dir}/${i}.sorted.rmdup.intergenic.unique.bam > ${sam_res_dir}/${i}.sorted.rmdup.intergenic.unique.paired_mapped.bam
#single read were mapped
#samtools view -b -F 4 -f 8 ${sam_res_dir}/${i}.sorted.rmdup.intergenic.unique.bam > ${sam_res_dir}/${i}.sorted.rmdup.intergenic.unique.one_read_mapped_1.bam
#samtools view -b -F 8 -f 4 ${sam_res_dir}/${i}.sorted.rmdup.intergenic.unique.bam > ${sam_res_dir}/${i}.sorted.rmdup.intergenic.unique.one_read_mapped_2.bam

#to obtain the unique reads mapped to intergenic region
samtools view -L ${intergenic} -h ${bam_dir}/${i}.sorted.bam| grep -P '^\@|NH:i:1\b'| samtools view -b > ${sam_res_dir}/${i}.sorted.intergenic.unique.bam

done

echo "gDNA_fragment_mapping.sh finished! "
