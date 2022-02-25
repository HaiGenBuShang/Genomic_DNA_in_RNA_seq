#!/bin/bash

data_dir=$(pwd)/../../data
file_list=${1}
#ref_dir=/path/to/HISAT2/index
ref_dir=/mnt/Miscrosoft/Wang_lab/Project/Reference/grch38
res_dir=$(pwd)/../../results/mapping
bam_out_dir=$(pwd)/../../results/mapping

#Trimmomatic file
#trimmo=/path/to/Trimmomatic/jar/file
trimmo=/home/hgbs/src/Trimmomatic-0.36/trimmomatic-0.36.jar
#adapter=/path/to/Trimmomatic/adapter/file
adapter=/home/hgbs/src/Trimmomatic-0.36/adapters/TruSeq3-PE.fa

#annotation file
#gencode_gtf=/path/to/gtf/file
gencode_gtf=/mnt/Miscrosoft/Shi_lab/genomicDNA/genomic_DNA_NB/gtf/gencode.v22.annotation.gtf
out_gtf_no_chr=$(pwd)/../../data/gencode.v22.annotation.no.chr.gtf
intergenic=$(pwd)/../../data/stringtie_merged_all_PA_and_RZ_0pct.intergenic.bed
anno=${out_gtf_no_chr}

./gDNA_mapping.sh ${file_list} ${data_dir} ${ref_dir} ${res_dir} ${trimmo} ${adapter} ${bam_out_dir}
./gDNA_stringtie.sh ${file_list} ${data_dir} ${ref_dir} ${res_dir} ${trimmo} ${adapter} ${anno} ${intergenic} ${bam_out_dir}

while [ ! -e "../intergenic_region/intergenic_region_info_wait_to_delete_automatically.txt" ]
do
sleep 2m
echo "Intergenic region bed file hasn't been generated! "
done
echo "Intergenic region has been generated! "

./gDNA_fragment_mapping.sh ${file_list} ${data_dir} ${ref_dir} ${res_dir} ${trimmo} ${adapter} ${anno} ${intergenic} ${bam_out_dir}
./gDNA_intergenic_ratio.sh ${file_list} ${data_dir} ${ref_dir} ${res_dir} ${trimmo} ${adapter} ${anno} ${intergenic} ${bam_out_dir}


echo "gDNA_estimate finished! " >> gDNA_finished_wait_to_delete_automatically.txt



