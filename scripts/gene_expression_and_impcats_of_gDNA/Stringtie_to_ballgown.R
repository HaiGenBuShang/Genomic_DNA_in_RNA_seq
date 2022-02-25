all_args <- commandArgs(trailingOnly = TRUE)

ballgown_dir <- all_args[[1]]

library(ballgown)
data_dir <- ballgown_dir
ballgown_expr <- ballgown(dataDir = data_dir, samplePattern = ".*",meas = "FPKM")
save(ballgown_expr,file="../../results/expression/ballgown_expr_data.Rdata")

no_merge_gene_fpkm <- gexpr(ballgown_expr)
save(no_merge_gene_fpkm,file="../../results/expression/no_merge_gene_fpkm.RData")

system("echo \"Stringtie_to_ballgown.R finished! \"")
