all_args <- commandArgs(trailingOnly = TRUE)

no_merge_gene_fpkm_RData <- all_args[[1]]


library(ballgown)
load(no_merge_gene_fpkm_RData)
no_merge_gene_fpkm <- log2(no_merge_gene_fpkm+0.01)

#results folder
constant <- 0.01
file_folder <- "../../results/expression/"

#calculate DEG
gene_expr <- no_merge_gene_fpkm
RZ_g_expr <- gene_expr[,31:48][,c(16:18,1:6,10:12,7:9,13:15)]
PA_pct_expr <- gene_expr[,1:30][,c(28:30,1:6,22:24,19:21,25:27)]

directory <- paste0(file_folder,"/DEGs_at_constant_",constant,"/")
dir.create(path = directory,recursive = TRUE,showWarnings = FALSE)

#RZ DEGs
ori_varianbles <- ls()

all_samples <- RZ_g_expr[,1:18]
figure_folder <- directory
file_name <- "RZ"
PA_and_RZ <- FALSE
source("./report_DEGs.R")
final_variables <- ls()
rm(list = c("final_variables",final_variables[!final_variables%in%ori_varianbles]))

message(paste0("RZ with contant: ",constant," finished!"))



#different genes
ori_varianbles <- ls()
all_samples <- PA_pct_expr[,1:18]
figure_folder <- directory
file_name <- "PA_pct"
PA_and_RZ <- FALSE
source("./report_DEGs.R")
final_variables <- ls()
rm(list = c("final_variables",final_variables[!final_variables%in%ori_varianbles]))

message(paste0("Poly (A) pct with contant: ",constant," finished!"))



#different genes
ori_varianbles <- ls()

all_samples <- cbind(PA_pct_expr[,1:18],RZ_g_expr[,1:18])
figure_folder <- directory
file_name <- "PA_VS_RZ"
PA_and_RZ <- TRUE
source("./report_DEGs.R")
final_variables <- ls()
rm(list = c("final_variables",final_variables[!final_variables%in%ori_varianbles]))

message(paste0("Poly (A) pct compared to Ribo-Zero with contant: ",constant," finished!"))


system("echo \"RZ_and_PA_DEGs.R finished! \"")