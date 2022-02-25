all_args <- commandArgs(trailingOnly = TRUE)

no_merge_gene_fpkm_RData <- all_args[[1]]

##function to select genes
select_genes <- function(expr_dt,expr_pct=0.7,expr_thres=0){
  n_clmn <- ncol(expr_dt)
  #if the expression of one gene in a samples is bigger than expr_thres
  total_count <- apply(expr_dt,1,function(x){
    sum(x > expr_thres)
  })
  total_count[rownames(expr_dt)] > expr_pct*n_clmn
}


############gene expression
library(ballgown)
load(no_merge_gene_fpkm_RData)

no_merge_gene_fpkm <- log2(no_merge_gene_fpkm+0.01)

##RZ gene expression
RZ_g_expr <- no_merge_gene_fpkm[,31:48]
RZ_g_expr <- RZ_g_expr[,c(16:18,1:6,10:12,7:9,13:15)]


kept_genes <- rownames(RZ_g_expr)[select_genes(expr_dt = RZ_g_expr,expr_pct = 0.3,
                                               expr_thres = log2(0.02))]



write.table(kept_genes,file = "../../results/expression/RZ_only_high_expressed_genes.txt",
            col.names = FALSE,row.names = FALSE,sep = "\t",quote = FALSE)

#PA pct expr
PA_pct_expr <- no_merge_gene_fpkm[,1:30][,c(28:30,1:6,22:24,19:21,25:27)]
PA_pct_kept_genes <- rownames(PA_pct_expr)[select_genes(expr_dt = PA_pct_expr,
                                                        expr_pct = 0.3,expr_thres = log2(0.02))]
write.table(PA_pct_kept_genes,file = "../../results/expression/PA_pct_only_high_expressed_genes.txt",
            col.names = FALSE,row.names = FALSE,sep = "\t",quote = FALSE)

#PA pct no and RZ no expr
PA_pct_no_VS_RZ_no_expr <- no_merge_gene_fpkm[,c(28:30,46:48)]
PA_pct_no_VS_RZ_no_kept_genes<-rownames(no_merge_gene_fpkm)[select_genes(expr_dt = PA_pct_no_VS_RZ_no_expr,
                                                                         expr_pct = 0.999, #in fact we kept genes,
                                                                         #which expression > log2(0.02)
                                                                         #in all no gDNA replicates
                                                                         expr_thres = log2(0.02))]
write.table(PA_pct_no_VS_RZ_no_kept_genes,
            file = "../../results/expression/PA_pct_no_VS_RZ_no_only_high_expressed_genes.txt",
            col.names = FALSE,row.names = FALSE,sep = "\t",quote = FALSE)

system("echo \"filter_low_expressed_genes.R finished! \"")