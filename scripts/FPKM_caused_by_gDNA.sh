#calculate FPKM caused by gDNA
cd residual_gDNA_and_expression_adjusting/
Rscript give_adjusting_expression_value_for_sequenced_samples.R "/path/to/mapping/directory/for/samples/to/be/gDNA/FPKM/estimated" "../../data/stringtie_merged_all_PA_and_RZ_0pct.intergenic.bed"
