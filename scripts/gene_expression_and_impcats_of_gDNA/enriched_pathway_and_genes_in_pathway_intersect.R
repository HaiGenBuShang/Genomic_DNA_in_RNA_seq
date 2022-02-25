all_args <- commandArgs(trailingOnly = TRUE)

PA_VS_RZ_RData <- all_args[[1]]
figure_dir <- all_args[[2]]


obtain_DEGs_compare <- function(DEG_RData,filter_gene_file){
  load(DEG_RData)
  filter_list <- read.table(filter_gene_file,
                            header = FALSE,sep = "\t",stringsAsFactors = FALSE)[,1]
  variables <- ls()
  all_DEGS <- lapply(all_DEGS,function(x){
    intersect(x,filter_list)
  })
  return(all_DEGS)
  rm(list = variables)
}


library(clusterProfiler)
warning("The varied version of clusterProfilter would change the number of enriched pathways!")

load(PA_VS_RZ_RData)

## all compared to PA with 0% gDNA
target_treatments <- c("Poly_A_No_VS_RZ_No","Poly_A_No_VS_RZ_0_01PCT","Poly_A_No_VS_RZ_0_1PCT",
                       "Poly_A_No_VS_RZ_1PCT","Poly_A_No_VS_RZ_10PCT")

target_treatments_DEGs <- all_sample_DEGs[target_treatments]

all_genes_considered <- read.table("../../results/expression/RZ_only_high_expressed_genes.txt",
                                   header = FALSE,sep = "\t",stringsAsFactors = FALSE)
PA_no_VS_RZ <- lapply(target_treatments_DEGs,function(x,y){
  intersect(rownames(x),y)
},y=all_genes_considered[,1])


# KEGG pathway enriching results
PA_no_VS_RZ_kegg <- lapply(PA_no_VS_RZ,function(x){
  g_id <- gsub("\\..*","",x)
  entrez_id <- bitr(geneID = g_id,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")$ENTREZID
  
  enrichKEGG(gene = entrez_id)
})

# genes mapped to enriched pathways
PA_no_VS_RZ_enriched_gene <- lapply(PA_no_VS_RZ_kegg,function(x){
  enrich_res <- DOSE:::get_enriched(x)@result
  if(nrow(enrich_res)>0){
    enriched_gene <- lapply(strsplit(enrich_res$geneID,split = "/"),function(x)x)
    names(enriched_gene) <- rownames(enrich_res)
  }else{
    enriched_gene <- list()
  }
  enriched_gene
})

# enriched pathways in each comparison combination
PA_no_VS_RZ_enriched_pathways <- lapply(PA_no_VS_RZ_enriched_gene,names)

VennDiagram::venn.diagram(x = PA_no_VS_RZ_enriched_pathways,
                          filename = paste0(figure_dir,"/PA_no_VS_RZ_kegg_venn_publish.tiff"),
                          fill=c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
                          category.names = c("0%","0.01%","0.1%","1%","10%"),
                          lwd=2,lty="blank",
                          cex=c(1,2,rep(1,6),2,rep(1,5),2,rep(1,8),2,rep(1,5),2,3),
                          fontface="bold",fontfamily="sans",
                          cat.fontfamily="sans",cat.fontface = "bold",cat.cex=2,
                          cat.pos=c(-50,-20,-112,112,20),cat.dist=c(0.2,0.2,0.22,0.2,0.2))






##Ribo-zero
RZ_0_01 <- obtain_DEGs_compare(DEG_RData = "../../results/expression/DEGs_at_constant_0.01/DEGs_RZ.RData",
                               filter_gene_file = "../../results/expression/RZ_only_high_expressed_genes.txt")


RZ_0_01_kegg <- lapply(RZ_0_01,function(x){
  g_id <- gsub("\\..*","",x)
  entrez_id <- bitr(geneID = g_id,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")$ENTREZID
  
  enrichKEGG(gene = entrez_id)
})


RZ_0_01_kegg_enriched_gene <- lapply(RZ_0_01_kegg,function(x){
  enrich_res <- DOSE:::get_enriched(x)@result
  if(nrow(enrich_res)>0){
    enriched_gene <- lapply(strsplit(enrich_res$geneID,split = "/"),function(x)x)
    names(enriched_gene) <- rownames(enrich_res)
  }else{
    enriched_gene <- list()
  }
  enriched_gene
})

##venn plot of PA_VS_RZ and RZ_VS_RZ enriched pathways
RZ_0_01_kegg_enriched_pathways <- lapply(RZ_0_01_kegg_enriched_gene,names)
PA_VS_RZ_and_RZ_enriched_pathways <- list(PA_VS_RZ=unique(do.call(c,PA_no_VS_RZ_enriched_pathways)),
                                          RZ=RZ_0_01_kegg_enriched_pathways$RZ_No_VS_RZ_10PCT)
VennDiagram::venn.diagram(x = PA_VS_RZ_and_RZ_enriched_pathways,
                          filename = paste0(figure_dir,"PA_VS_RZ_and_RZ_kegg_venn_publish.tiff"),
                          fill=c("cornflowerblue", "darkorchid1"),
                          category.names = c("PA VS RZ","RZ"),
                          lwd=2,lty="blank",
                          cex=2,cat.pos=c(190,170),
                          fontface="bold",fontfamily="sans",
                          cat.fontfamily="sans",cat.fontface = "bold",cat.cex=2,direct.area = FALSE,ext.text=FALSE)


## most enriched pathway
the_pathway <- RZ_0_01_kegg$RZ_No_VS_RZ_10PCT@result[1,1]
pathway_enriched_genes <- list(PA_VS_RZ=sapply(strsplit(PA_no_VS_RZ_kegg$Poly_A_No_VS_RZ_10PCT@result[the_pathway,"geneID"],"/"),
                                               function(x)x)[,1],
                               RZ=sapply(strsplit(RZ_0_01_kegg$RZ_No_VS_RZ_10PCT@result[the_pathway,"geneID"],"/"),
                                         function(x)x)[,1])

VennDiagram::venn.diagram(x = pathway_enriched_genes,
                          filename = paste0(figure_dir,
                                            "PA_VS_RZ_and_RZ_kegg_",the_pathway,"_venn_publish.tiff"),
                          fill=c("cornflowerblue", "darkorchid1"),
                          category.names = c("PA VS RZ","RZ"),
                          lwd=2,lty="blank",
                          cex=2,cat.pos=c(-40,50),cat.dist=c(0.04,0.025),
                          fontface="bold",fontfamily="sans",
                          cat.fontfamily="sans",cat.fontface = "bold",cat.cex=2,direct.area = FALSE,ext.text=FALSE)

RZ_0_01_kegg$RZ_No_VS_RZ_10PCT@result[1,3:4]
PA_no_VS_RZ_kegg$Poly_A_No_VS_RZ_10PCT@result[the_pathway,3:4]

pathway_matrix <- apply(RZ_0_01_kegg$RZ_No_VS_RZ_10PCT@result[1,3:4],2,function(x)as.numeric(strsplit(x,"/")[[1]]))
RZ_0_01_the_pathway_matrix <- matrix(c(pathway_matrix[1,1],pathway_matrix[1,2]-pathway_matrix[1,1],
                                       pathway_matrix[2,1]-pathway_matrix[1,1],
                                       pathway_matrix[2,2]-pathway_matrix[1,2]-pathway_matrix[2,1]+pathway_matrix[1,1]),
                                     nrow = 2,byrow = TRUE)
dimnames(RZ_0_01_the_pathway_matrix) <- list(c("In_Pathway","Not_in_Pathway"),c("DEG","Non-DEG"))

write.table(addmargins(RZ_0_01_the_pathway_matrix),
            file = paste0("../../results/expression/RZ_pathway_",the_pathway,"_enriched_results.txt"),
            col.names = TRUE,row.names = TRUE,sep = "\t",quote = FALSE)

pathway_matrix <- apply(PA_no_VS_RZ_kegg$Poly_A_No_VS_RZ_10PCT@result[the_pathway,3:4],2,
                        function(x)as.numeric(strsplit(x,"/")[[1]]))
PA_no_vs_RZ_pathway_matrix<-matrix(c(pathway_matrix[1,1],pathway_matrix[1,2]-pathway_matrix[1,1],
                                     pathway_matrix[2,1]-pathway_matrix[1,1],
                                     pathway_matrix[2,2]-pathway_matrix[1,2]-pathway_matrix[2,1]+pathway_matrix[1,1]),
                                   nrow = 2,byrow = TRUE)
dimnames(PA_no_vs_RZ_pathway_matrix) <- list(c("In_Pathway","Not_in_Pathway"),c("DEG","Non-DEG"))
write.table(addmargins(PA_no_vs_RZ_pathway_matrix),
            file = paste0("../../results/expression/PA_no_vs_RZ_pathway_",the_pathway,"_enriched_results.txt"),
            col.names = TRUE,row.names = TRUE,sep = "\t",quote = FALSE)



a <- PA_no_VS_RZ_kegg$Poly_A_No_VS_RZ_10PCT@result$pvalue
all_adjusted_p<- c()
for(i in sum(PA_no_vs_RZ_pathway_matrix[,1]):(PA_no_vs_RZ_pathway_matrix[1,1]+1)){
  total_DEGs_in_bg <- i
  ov_gene <- PA_no_vs_RZ_pathway_matrix[1,1]
  
  the_pathway_index <- grep(the_pathway,rownames(PA_no_VS_RZ_kegg$Poly_A_No_VS_RZ_10PCT@result))
  the_pathway_value <- PA_no_VS_RZ_kegg$Poly_A_No_VS_RZ_10PCT@result[the_pathway_index,c("GeneRatio","BgRatio")]
  the_pathway_value <- lapply(apply(the_pathway_value[1,],1,strsplit,split="/"),function(x)lapply(x,as.integer))
  
  
  
  p_v <- fisher.test(matrix(c(ov_gene,total_DEGs_in_bg-ov_gene,
                              the_pathway_value[[1]][[2]][1]-ov_gene,
                              the_pathway_value[[1]][[2]][2]-total_DEGs_in_bg-the_pathway_value[[1]][[2]][1]+ov_gene),
                            nrow = 2,byrow = FALSE),
                     alternative = "greater")$p.value
  
  a[the_pathway_index] <- p_v
  # message(paste0("total_DEGs_in_bg = ",i))
  # print(p.adjust(a,method = "bonferroni")[the_pathway_index])
  adjusted_p <- p.adjust(a,method = "BH")
  all_adjusted_p <- c(all_adjusted_p,adjusted_p[the_pathway_index])
}
names(all_adjusted_p) <- sum(PA_no_vs_RZ_pathway_matrix[,1]):(PA_no_vs_RZ_pathway_matrix[1,1]+1)


pdf(paste0(figure_dir,"/P_value_and_DEG_number.pdf"))
par(mai=par("mai")+0.2)
plot(all_adjusted_p,xlab="DEGs in Background",ylab="Adjusted P-Value",xaxt="n",cex.lab=1.5,font.sub=2,cex.axis = 1.5,
     main="Adjusted P-value at\nDifferent DEGs in Background",cex.main=1.8)
axis(side = 1,
     at = floor(quantile((PA_no_vs_RZ_pathway_matrix[1,1]+1):sum(PA_no_vs_RZ_pathway_matrix[,1]))-PA_no_vs_RZ_pathway_matrix[1,1]),
     labels = (sum(PA_no_vs_RZ_pathway_matrix[,1]):(PA_no_vs_RZ_pathway_matrix[1,1]+1))[
       floor(quantile((PA_no_vs_RZ_pathway_matrix[1,1]+1):sum(PA_no_vs_RZ_pathway_matrix[,1]))-PA_no_vs_RZ_pathway_matrix[1,1])],cex.axis=1.5)
abline(h = 0.05,lty=2)
box()
dev.off()


# DEG numbers between PA and RZ
pdf(paste0(figure_dir,"/number_of_DEGs_RZ_vs_PA_and_RZ_vs_RZ.pdf"),
    width = 6.5,height = 6.5, colormodel = "cmyk")
DEG_matrix <- c(sapply(PA_no_VS_RZ,length),length(RZ_0_01[[length(RZ_0_01)]]))
coor <- barplot(as.vector(DEG_matrix),col = c("#377EB8","#377EB8","#377EB8","#377EB8","#377EB8","#E41A1C"),
                main = "DEGs between Ribo-Zero and\nPoly (A) Selection/Ribo-Zero",
                cex.main=1.8,sub = "gDNA Concentration",
                cex.sub=1.8,font.sub=2,cex.axis = 1.5,ylim = c(0,max(DEG_matrix)+1000),
                border = c("#377EB8","#377EB8","#377EB8","#377EB8","#377EB8","#E41A1C"),yaxt="n")
axis(2,at = quantile(0:floor((max(DEG_matrix)+1000)/1000),probs=seq(0,1,0.25))*1000,cex.axis=1.5)

mtext(text = c("0%","0.01%","0.1%","1%","10%","10%"),side = 1,line = 1.2,
      at = coor[,1],font = 2,cex = 1.4,padj = 0.5)

legend("topleft",legend = c("RZ vs PA with 0% gDNA","RZ vs RZ with 0% gDNA"),
       fill = c("#377EB8","#E41A1C"),bty = "n",cex = 1.5,border = c("#377EB8","#E41A1C"))
dev.off()



PA_no_VS_RZ_kegg_enriched <- lapply(PA_no_VS_RZ_kegg,DOSE:::get_enriched)

PA_no_VS_RZ_kegg_gene_in_background <- sapply(PA_no_VS_RZ_kegg_enriched,function(x){
  unique(as.numeric(gsub(".*/","",x@result$GeneRatio)))
})

# DEG numbers in background
pdf(paste0(figure_dir,"/number_of_DEGs_in_bg_RZ_vs_PA_and_RZ_vs_RZ.pdf"),
    width = 6.5,height = 6.5, colormodel = "cmyk")
DEG_matrix <- c(PA_no_VS_RZ_kegg_gene_in_background,
                as.numeric(strsplit(RZ_0_01_kegg$RZ_No_VS_RZ_10PCT@result[1,]$GeneRatio,split = "/")[[1]][2]))
coor <- barplot(as.vector(DEG_matrix),col = c("#377EB8","#377EB8","#377EB8","#377EB8","#377EB8","#E41A1C"),
                main = "DEGs in the Background of\nEnrichment Analysis",cex.main=1.8,
                sub = "gDNA Concentration",
                cex.sub=1.8,font.sub=2,cex.axis = 1.5,ylim = c(0,max(as.vector(DEG_matrix))+1000),
                border = c("#377EB8","#377EB8","#377EB8","#377EB8","#377EB8","#E41A1C"),yaxt="n")
axis(2,at = quantile(0:floor((max(DEG_matrix)+1000)/1000),probs=seq(0,1,0.25))*1000,cex.axis=1.5)
mtext(text = c("0%","0.01%","0.1%","1%","10%","10%"),side = 1,line = 1.2,
      at = coor[,1],font = 2,cex = 1.5,padj = 0.5)
legend("topleft",legend = c("RZ vs PA with 0% gDNA","RZ vs RZ with 0% gDNA"),
       fill = c("#377EB8","#E41A1C"),bty = "n",cex = 1.5,border = c("#377EB8","#E41A1C"))
dev.off()

system("echo \"enriched_pathway_and_genes_in_pathway_intersect.R finished! \"")