args <- commandArgs(trailingOnly = TRUE)

mapping_dir <- args[[1]]
intergenic_bed_file <- args[[2]]
  

get_gDNA_FPKM <- function(mapped_file,target_mapped_file,intergenic_bed_file){
  cat("\033[48;5;9;1mThis function is only suitable for Ribo-Zero!\033[0m\n")
  
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
setwd(mapping_dir)
total_mapped_file <- list.files(path = "./",pattern = "sorted.bam.stats.txt",recursive = TRUE)
intergenic_mapped_file <- list.files(path = "./",pattern = "sorted.intergenic.unique.bam.stats.txt",recursive = TRUE)
all_sample_gDNA_FPKM <- mapply(get_gDNA_FPKM,mapped_file=total_mapped_file,target_mapped_file=intergenic_mapped_file,
                               intergenic_bed_file=intergenic_bed_file)
write.table(all_sample_gDNA_FPKM,
            file = paste0(current_dir,"/../../results/expression/gDNA_FPKM_for_sequenced_samples.txt"),
            col.names = FALSE,row.names = names(all_sample_gDNA_FPKM),sep = "\t",quote = FALSE)
