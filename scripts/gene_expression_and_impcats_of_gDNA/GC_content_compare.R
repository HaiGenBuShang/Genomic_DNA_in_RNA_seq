all_args <- commandArgs(trailingOnly = TRUE)

RZ_cor_gene_file <- all_args[[1]]
RZ_only_high_expressed_gene_file <- all_args[[2]]

library(biomaRt)
RZ_high_expressed_genes <- read.table(RZ_only_high_expressed_gene_file,
                                      header = FALSE,sep = "\t",stringsAsFactors = FALSE)[,1]
RZ_cor_genes <- read.table(RZ_cor_gene_file,
                           header = FALSE,sep = "\t",stringsAsFactors = FALSE)[,1]


mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",host = "http://uswest.ensembl.org")

all_gene_GC_content <- getBM(attributes = c("percentage_gene_gc_content","ensembl_gene_id_version"),
                             filters = "ensembl_gene_id_version",
                             values = RZ_high_expressed_genes,
                             mart = mart,useCache = FALSE)

rownames(all_gene_GC_content) <- all_gene_GC_content[,"ensembl_gene_id_version"]

cor_gene_gc_content <- all_gene_GC_content[RZ_cor_genes,1]
other_gene_gc_content <- all_gene_GC_content[setdiff(all_gene_GC_content$ensembl_gene_id_version,RZ_cor_genes),1]

message(paste0("The number of genes correlated with gDNA have a GC content value: ",sum(!is.na(cor_gene_gc_content)),
               "\nThe number of genes not correlated with gDNA have a GC content value: ",length(other_gene_gc_content),
               "\nOnly high expressed genes were counted!"))
t.test(cor_gene_gc_content,other_gene_gc_content)

system("echo \"GC_content_compare.R finished! \" ")
