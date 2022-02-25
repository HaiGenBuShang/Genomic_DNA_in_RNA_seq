#check if all gDNA estimated program finished
gDNA_estimate_numbers=$(grep "gDNA_estimate.sh" ../Reads_mapping_and_gDNA_estimate.sh|wc -l)
gDNA_estimate_finished_file=gDNA_finished_wait_to_delete_automatically.txt
finished_gDNA=0

#sleep 4h

while [ "${finished_gDNA}" -lt "${gDNA_estimate_numbers}" ]
do
sleep 1m
if [ -e "${gDNA_estimate_finished_file}" ] ; then
finished_gDNA=$(wc -l < ${gDNA_estimate_finished_file})
fi
printf "gDNA_estimate.sh finished number: ${finished_gDNA}!\n"
done
echo "All gDNA_estimate.sh finished!"
rm ${gDNA_estimate_finished_file} ../intergenic_region/intergenic_region_info_wait_to_delete_automatically.txt
