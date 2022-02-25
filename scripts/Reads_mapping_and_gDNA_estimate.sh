#!/bin/bash

cd intergenic_region/
./Run_produce_no.chr.gtf.sh > no_chr_gtf.log 2>&1 &
cd - > /dev/null

cd mapping_to_genome_and_transcripts_assembling/
nohup ./gDNA_estimate.sh gDNA_file.list00 > gDNA_estimate_0.log 2>&1 &
nohup ./gDNA_estimate.sh gDNA_file.list01 > gDNA_estimate_1.log 2>&1 &
nohup ./gDNA_estimate.sh gDNA_file.list02 > gDNA_estimate_2.log 2>&1 &
nohup ./gDNA_estimate.sh gDNA_file.list03 > gDNA_estimate_3.log 2>&1 &
nohup ./gDNA_estimate_finish.sh > gDNA_finished_or_not.log 2>&1 &
cd - > /dev/null

# file_to_merge_gtf.list
cd intergenic_region/
nohup ./produce_intergenic_region_bed.sh > gDNA_intergenic.log 2>&1 &
cd - > /dev/null

