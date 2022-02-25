#!/bin/bash

gencode_gtf=/mnt/Miscrosoft/Shi_lab/genomicDNA/genomic_DNA_NB/gtf/gencode.v22.annotation.gtf
out_gtf_no_chr=$(pwd)/../../data/gencode.v22.annotation.no.chr.gtf
out_gene_type_and_name=$(pwd)/../../data/gene_type_and_gene_name.txt
Rscript produce_no.chr.gtf.R ${gencode_gtf} ${out_gtf_no_chr}

awk -F ';' '{print $1,$2,$4}' ${gencode_gtf} |grep type|awk -F '\t' '{print $9,$10,$11}'|awk -F ' ' '{print $2"\t"$4"\t"$6}' > ${out_gene_type_and_name}
