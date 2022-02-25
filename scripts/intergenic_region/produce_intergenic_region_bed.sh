#!/bin/bash

ref_gtf=../../data/gencode.v22.annotation.no.chr.gtf
out_gtf_file=../../data/stringtie_all_PA_and_RZ_0pct_merged.gtf
out_bed_file=../../data/stringtie_merged_all_PA_and_RZ_0pct.intergenic.bed
gtf_list=all_PA_pct_and_RZ_0pct_gtf.list
bam_out_dir=$(pwd)/../../results/mapping

bam=$(find ${bam_out_dir} -name "*bam"|head -n 1)
merged_gtf=${out_gtf_file}




#check if to start intergenic region produce
stringtie_numbers=$(grep "gDNA_estimate.sh" ../Reads_mapping_and_gDNA_estimate.sh|wc -l)
stringtie_file=../mapping_to_genome_and_transcripts_assembling/wait_to_delete_automatically.txt
finished_stringtie=0

#sleep 5h

while [ "${finished_stringtie}" -lt "${stringtie_numbers}" ]
do
sleep 1m
if [ -e "${stringtie_file}" ] ; then
finished_stringtie=$(wc -l < ../mapping_to_genome_and_transcripts_assembling/wait_to_delete_automatically.txt)
fi
printf "waiting for StringTie finishing! \n"
echo $finished_stringtie
done
rm ../mapping_to_genome_and_transcripts_assembling/wait_to_delete_automatically.txt
echo "finised $finished_stringtie. Start merging! "
sleep 1m

##******************
##To merge gtf files
##******************
stringtie --merge -G ${ref_gtf} -o ${out_gtf_file} ${gtf_list}


##********************************
##To produce intergenic region bed
##********************************
#to produce the genome length file for bed complement
samtools view -H ${bam} > tmp.chr.length.txt
#echo $bam
grep SQ  tmp.chr.length.txt | cut -d ':' -f 2,3 | sed -s 's/LN://g' > gencode.v22.anno.no.chr.length.txt
rm tmp.chr.length.txt

awk -F '\t' '($3=="transcript") {printf("%s\t%d\t%s\n",$1,int($4)-1,$5);}' ${merged_gtf} |sort -t $'\t' -k1,1 -k2,2n| bedtools merge > stringtie_merged_all_PA_and_RZ_0pct.no.chr.bed
bedtools complement -i stringtie_merged_all_PA_and_RZ_0pct.no.chr.bed -g gencode.v22.anno.no.chr.length.txt > ${out_bed_file}
rm gencode.v22.anno.no.chr.length.txt stringtie_merged_all_PA_and_RZ_0pct.no.chr.bed

printf "Intergenic_bed finished! \n"
echo "intergenic_region_bed_finished! " > intergenic_region_info_wait_to_delete_automatically.txt
