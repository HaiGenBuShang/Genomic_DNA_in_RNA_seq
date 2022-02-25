args <- commandArgs(trailingOnly = TRUE)

gDNA_project_mapping_dir <- args[[1]]
# TNBC_sample_mapping_dir <- args[[2]]
# Scientific_data_mapping_dir <- args[[3]]
# intergenic_bed_file <- args[[4]]
intergenic_bed_file <- args[[2]]
  

get_gDNA_FPKM <- function(mapped_file,target_mapped_file,intergenic_bed_file){
  message("\033[48;5;9;1mThis function is only suitable for Ribo-Zero!\033[0m")
  
  system(command = paste0("cat ",intergenic_bed_file," | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' > tmp_total_intergenic_length.txt"))
  intergenic_length <- read.table("tmp_total_intergenic_length.txt",header = FALSE,stringsAsFactors = FALSE)[1,1]
  system(command = "rm tmp_total_intergenic_length.txt")
  
  target_mapped_reads <- obtain_mapped_reads(target_mapped_file)
  target_mapped_paired_reads <- obtain_mapped_paired_reads(target_mapped_file)
  
  total_reads <- obtain_mapped_reads(mapped_file)
  

  target_one_read_mapped <- target_mapped_reads - target_mapped_paired_reads
  gDNA_intergenic_frament <- target_one_read_mapped + target_mapped_paired_reads/2
  (FPKM_gDNA <- (gDNA_intergenic_frament/(total_reads*intergenic_length))*10e9)
  
}


obtain_mapped_reads <- function(mapped_file){
  reads_mapping_info <- read.table(mapped_file,header = FALSE,sep = "\n",stringsAsFactors = FALSE)
  mapped_reads <- grep("reads mapped:",reads_mapping_info[,1],value = TRUE)
  total_mapped_reads <- strsplit(mapped_reads,split = "\t")[[1]][3]
  as.numeric(total_mapped_reads)
}

obtain_mapped_paired_reads <- function(mapped_file){
  reads_mapping_info <- read.table(mapped_file,header = FALSE,sep = "\n",stringsAsFactors = FALSE)
  mapped_reads <- grep("reads mapped and paired:",reads_mapping_info[,1],value = TRUE)
  total_mapped_reads <- strsplit(mapped_reads,split = "\t")[[1]][3]
  as.numeric(total_mapped_reads)
}

target_mapping_ratio <- function(total_mapped_file,target_mapped_file){
  total_mapped_reads <- obtain_mapped_reads(total_mapped_file)
  target_mapped_reads <- obtain_mapped_reads(target_mapped_file)
  target_mapped_reads/total_mapped_reads
}


current_dir <- getwd()
setwd(gDNA_project_mapping_dir)
total_mapped_file <- list.files(path = "./",pattern = "sorted.bam.stats.txt",recursive = TRUE)
intergenic_mapped_file <- list.files(path = "./",pattern = "sorted.intergenic.unique.bam.stats.txt",recursive = TRUE)

##RZ gDNA FPKM
RZ_total_mapped_file <- grep("RZ",total_mapped_file,value = TRUE)
RZ_intergenic_mapped_file <- grep("RZ",intergenic_mapped_file,value = TRUE)

all_RZ_gDNA_FPKM <- mapply(get_gDNA_FPKM,mapped_file=RZ_total_mapped_file,target_mapped_file=RZ_intergenic_mapped_file,
                           intergenic_bed_file=intergenic_bed_file)

names(all_RZ_gDNA_FPKM) <- paste0(rep(c("0.01_pct","0.1_pct","10_pct","1_pct","No_DNase","0_pct"),each=3),"_",1:3)
write.table(all_RZ_gDNA_FPKM,
            file = paste0(current_dir,"/../../results/expression/RZ_gDNA_FPKM_value.txt"),
            col.names = FALSE,row.names = names(all_RZ_gDNA_FPKM),sep = "\t",quote = FALSE)


##PA gDNA FPKM
PA_total_mapped_file <- grep("PCT.*_[123].*/|Native.*_[123]|No.*_[123]",total_mapped_file,value = TRUE)
PA_intergenic_mapped_file <- grep("PCT.*_[123].*/|Native.*_[123]|No.*_[123]",intergenic_mapped_file,value = TRUE)

all_PA_gDNA_FPKM <- mapply(get_gDNA_FPKM,mapped_file=PA_total_mapped_file,target_mapped_file=PA_intergenic_mapped_file,
                           intergenic_bed_file=intergenic_bed_file)
names(all_PA_gDNA_FPKM) <- paste0(rep(c("0.01_pct","0.1_pct","10_pct","1_pct","No_DNase","0_pct"),each=3),"_",1:3)
write.table(all_PA_gDNA_FPKM,file = paste0(current_dir,"/../../results/expression/PA_gDNA_FPKM_value.txt"),
            col.names = FALSE,row.names = names(all_PA_gDNA_FPKM),sep = "\t",quote = FALSE)

# ##TNBC gDNA FPKM
# setwd(TNBC_sample_mapping_dir)
# TNBC_mapped_file <- list.files(path = "./",pattern = "sorted.rmdup.bam.stats.txt",recursive = TRUE)
# TNBC_interg_mapped_file <- list.files(path = "./",pattern = "sorted.rmdup.intergenic.unique.bam.stats.txt",
#                                       recursive = TRUE)
# all_TNBC_gDNA_FPKM <- mapply(get_gDNA_FPKM,mapped_file=TNBC_mapped_file,target_mapped_file=TNBC_interg_mapped_file,
#                              intergenic_bed_file=intergenic_bed_file)
# message("\n\n\033[48;5;9;1mTNBC gDNA FPKM\033[0m\n")
# all_TNBC_gDNA_FPKM
# 
# ##Scientific Data gDNA FPKM
# setwd(Scientific_data_mapping_dir)
# SD_mapped_file <- list.files(path = "./",pattern = "sorted.rmdup.bam.stats.txt",recursive = TRUE)
# SD_interg_mapped_file <- list.files(path = "./",pattern = "sorted.rmdup.intergenic.unique.bam.stats.txt",recursive = TRUE)
# all_SD_gDNA_FPKM <- mapply(get_gDNA_FPKM,mapped_file=SD_mapped_file,target_mapped_file=SD_interg_mapped_file,
#                              intergenic_bed_file=intergenic_bed_file)
# message("\n\n\033[48;5;9;1mScientific Data gDNA FPKM\033[0m\n")
# all_SD_gDNA_FPKM

system("echo \"give_adjusting_expression_value.R\" finished! ")