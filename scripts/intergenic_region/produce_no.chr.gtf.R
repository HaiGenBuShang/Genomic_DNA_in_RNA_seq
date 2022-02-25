args <- commandArgs(trailingOnly = TRUE)

gtf_file_with_chr <- args[[1]]
out_gtf_without_chr <- args[[2]]

if(!file.exists("../../data/gencode.v22.annotation.no.chr.gtf")){
  gtf_with_chr <- read.table(gtf_file_with_chr,header = FALSE,sep = "\t",stringsAsFactors = FALSE,quote = "")
  gtf_with_chr[,1] <- gsub("chr","",gtf_with_chr[,1])
  gtf_with_chr[,1] <- gsub("M","MT",gtf_with_chr[,1])
  
  gtf_header <- read.table(gtf_file_with_chr,header = FALSE,sep = "\t",stringsAsFactors = FALSE,nrows = 5,comment.char = "")
  gtf_header[6,1] <- "##drop \"chr\" by Xiangnan and change \"chrM\" to \"MT\", on 20210912"
  write.table(gtf_header,file = out_gtf_without_chr,sep = "\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
  write.table(gtf_with_chr,file = "tmp_gencode.v22.gtf",sep = "\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
  
  system(paste0("cat tmp_gencode.v22.gtf >> ",out_gtf_without_chr))
  system("rm tmp_gencode.v22.gtf")
  system("echo \"no.chr.gtf finished! \" ")
}
