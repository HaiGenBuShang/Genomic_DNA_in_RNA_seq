args <- commandArgs(trailingOnly = TRUE)

gDNA_project_mapping_dir <- args[[1]]
gDNA_to_be_estimate_dir <- args[[2]]



#To infer residual gDNA
current_dir <- getwd()
setwd(gDNA_project_mapping_dir)


total_mapped_file <- list.files(path = "./",pattern = "sorted.bam.stats.txt",recursive = TRUE)
intergenic_mapped_file <- list.files(path = "./",pattern = "sorted.intergenic.unique.bam.stats.txt",recursive = TRUE)

obtain_mapped_reads <- function(mapped_file){
  reads_mapping_info <- read.table(mapped_file,header = FALSE,sep = "\n",stringsAsFactors = FALSE)
  mapped_reads <- grep("reads mapped:",reads_mapping_info[,1],value = TRUE)
  total_mapped_reads <- strsplit(mapped_reads,split = "\t")[[1]][3]
  as.numeric(total_mapped_reads)
}

target_mapping_ratio <- function(total_mapped_file,target_mapped_file){
  total_mapped_reads <- obtain_mapped_reads(total_mapped_file)
  target_mapped_reads <- obtain_mapped_reads(target_mapped_file)
  target_mapped_reads/total_mapped_reads
}

if(length(total_mapped_file)==0|length(intergenic_mapped_file)==0){
  all_mapping_ratio <- read.table(paste0(current_dir,"../../data/all_mapping_ratio.txt"),
                                  header = FALSE,sep = "\t",stringsAsFactors = FALSE,row.names = 1)[,1]
}else{
  all_mapping_ratio <- mapply(target_mapping_ratio,total_mapped_file=total_mapped_file,
                              target_mapped_file=intergenic_mapped_file)
  write.table(all_mapping_ratio,file = paste0(current_dir,"../../data/all_mapping_ratio.txt"),
              col.names = FALSE,row.names = names(all_mapping_ratio),quote = FALSE,sep = "\t")
}


mapping_rate_SE <- all_mapping_ratio




gDNA_concentration <- c(rep(c(0,0.0001,0.001,0.01,0.1),each=3))
gDNA_concentration <- gDNA_concentration/(1+gDNA_concentration)

mapping_rate_SE_PA_pct <- mapping_rate_SE[c(43:45,1:3,7:9,31:33,25:27)]
PA_fit_SE <- lm(mapping_rate_SE_PA_pct~gDNA_concentration)
message("\n\n\033[48;5;9;1mPoly (A) Selection linear regresion results\033[0m\n")
summary(PA_fit_SE)
library(lm.beta)
lm.beta(PA_fit_SE)

mapping_rate_SE_RZ_pct <- mapping_rate_SE[c(43:45,1:3,7:9,31:33,25:27)+3]
RZ_fit_SE <- lm(mapping_rate_SE_RZ_pct~gDNA_concentration)
message("\n\n\033[48;5;9;1mRibo-Zero linear regresion results\033[0m\n")
summary(RZ_fit_SE)
lm.beta(RZ_fit_SE)




#all annotated genes
gencode_v22_anno <- read.table("../../data/gene_type_and_gene_name.txt",header = FALSE,
                               stringsAsFactors = FALSE,sep = "\t")
total_annotated_genes <- nrow(gencode_v22_anno)
protein_coding_annotated <- table(gencode_v22_anno[,2])["protein_coding"]
non_protein_coding_annotated <- total_annotated_genes-protein_coding_annotated


#residual gDNA
residual_gDNA_in_total_RNA_SE <- (RZ_fit_SE$coefficients[1] - PA_fit_SE$coefficients[1]-
                                    PA_fit_SE$coefficients[1]*(non_protein_coding_annotated/protein_coding_annotated)
)/RZ_fit_SE$coefficients[2]
RZ_unannotated_time_alpha <- (PA_fit_SE$coefficients[1]+
                                PA_fit_SE$coefficients[1]*(non_protein_coding_annotated/protein_coding_annotated))

print(paste0("residual gDNA in total RNA: ",residual_gDNA_in_total_RNA_SE))
print(paste0("alpha * cDNAIR_RZ: ", RZ_unannotated_time_alpha))



#estimating gDNA contamination
setwd(gDNA_to_be_estimate_dir)
TNBC_mapped_file <- list.files(path = "./",pattern = "sorted.bam.stats.txt",recursive = TRUE)
TNBC_interg_mapped_file <- list.files(path = "./",pattern = "sorted.intergenic.unique.bam.stats.txt",recursive = TRUE)

TNBC_mapping_ratio <- mapply(target_mapping_ratio,total_mapped_file=TNBC_mapped_file,
                             target_mapped_file=TNBC_interg_mapped_file)

#names(TNBC_mapping_ratio) <- NULL
print(sapply(TNBC_mapping_ratio,function(x){
  (x - RZ_unannotated_time_alpha)/RZ_fit_SE$coefficients[2]
}))


# setwd("/mnt/Miscrosoft/Shi_lab/genomicDNA/genomic_DNA_NB/Scientific_data_mapping/merge_gtf_results/")
# SD_mapped_file <- list.files(path = "./",pattern = "sorted.bam.stats.txt",recursive = TRUE)
# SD_interg_mapped_file <- list.files(path = "./",pattern = "sorted.intergenic.unique.bam.stats.txt",recursive = TRUE)
# 
# SD_mapping_ratio <- mapply(target_mapping_ratio,total_mapped_file=SD_mapped_file,
#                              target_mapped_file=SD_interg_mapped_file)
# #names(SD_mapping_ratio) <- NULL
# sapply(SD_mapping_ratio,function(x){
#   (x - RZ_unannotated_time_alpha)/RZ_fit_SE$coefficients[2]
# })


# write.table(round(sapply(SD_mapping_ratio,function(x){
#   (x - RZ_unannotated_time_alpha)/RZ_fit_SE$coefficients[2]
# })*100,digits = 2),file = "~/gDNA.txt",quote = FALSE,col.names = FALSE,row.names = FALSE)


