printf "\033[35;1mYou may first run \"Reads_mapping_and_gDNA_estimate.sh\" first! \033[0m\n"

cd $(pwd)/gene_expression_and_impcats_of_gDNA
#./mk_ln_for_ballgown.sh "../../results/mapping/BallgownTable/" "../../results/ballgown_for_R/" "all_ballgown_dir.txt"

#Rscript ./install_R_packages.R

#Rscript ./Stringtie_to_ballgown.R "../../results/ballgown_for_R/"

#Rscript ./filter_low_expressed_genes.R "../../results/expression/no_merge_gene_fpkm.RData"

#Rscript ./PCA_and_heatmap.R "../../results/expression/no_merge_gene_fpkm.RData" "../../results/expression/figure/"

#Rscript ./RZ_and_PA_DEGs.R "../../results/expression/no_merge_gene_fpkm.RData"

#Rscript ./cor_genes_and_DEG_number_kegg_and_density.R "../../results/expression/no_merge_gene_fpkm.RData" "../../results/expression/figure/"

#Rscript ./enriched_pathway_and_genes_in_pathway_intersect.R "../../results/expression/DEGs_at_constant_0.01/DEGs_PA_VS_RZ.RData" "../../results/expression/figure/"

cd - > /dev/null

cd $(pwd)/residual_gDNA_and_expression_adjusting
Rscript ./give_adjusting_expression_value.R "../../results/mapping/" "../../data/stringtie_merged_all_PA_and_RZ_0pct.intergenic.bed"

cd - > /dev/null

cd $(pwd)/gene_expression_and_impcats_of_gDNA
Rscript ./Adjusting_FPKM_by_substract_gDNA_FPKM.R "../../results/expression/no_merge_gene_fpkm.RData" "../../results/expression/figure/"

#Rscript ./GC_content_compare.R "../../results/expression/RZ_genes_correlated_with_gDNA.txt" "../../results/expression/RZ_only_high_expressed_genes.txt"
cd - > /dev/null
