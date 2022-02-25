all_args <- commandArgs(trailingOnly = TRUE)

no_merge_gene_fpkm_RData <- all_args[[1]]
figure_dir <- all_args[[2]]


##obtain DEGs between each replicates
obtain_DEGs_compare <- function(DEG_RData,filter_gene_file){
  load(DEG_RData)
  filter_list <- read.table(filter_gene_file,
                            header = FALSE,sep = "\t",stringsAsFactors = FALSE)[,1]
  variables <- ls()
  all_DEGS <- lapply(all_DEGS,function(x){
    intersect(x,filter_list)
  })
  #return(sapply(all_DEGS,length))
  return(all_DEGS)
  rm(list = variables)
}

#function to help draw density plot
compare_two_density_y <- function(density_list){
  if(length(density_list)==1){
    return(max(density_list[[1]]$y))
  }else{
    max_1 <- max(density_list[[1]]$y)
    max(c(max_1,compare_two_density_y(density_list[-1])))
  }
}

compare_two_density_x <- function(density_list){
  if(length(density_list)==1){
    return(max(density_list[[1]]$x))
  }else{
    max_1 <- max(density_list[[1]]$x)
    max(c(max_1,compare_two_density_x(density_list[-1])))
  }
}


compare_two_density_x_2 <- function(density_list){
  if(length(density_list)==1){
    return(min(density_list[[1]]$x))
  }else{
    min_1 <- min(density_list[[1]]$x)
    min(c(min_1,compare_two_density_x(density_list[-1])))
  }
}


library(ballgown)
library(clusterProfiler)
load(no_merge_gene_fpkm_RData)
no_merge_gene_fpkm <- log2(no_merge_gene_fpkm+0.01)

RZ_g_expr <- no_merge_gene_fpkm[,31:48][,c(16:18,1:6,10:12,7:9,13:15)]
PA_pct_expr <- no_merge_gene_fpkm[,1:30][,c(28:30,1:6,22:24,19:21,25:27)]

#Ribo-Zero genes correlated with gDNA
RZ_expred_genes <- read.table(file = "../../results/expression/RZ_only_high_expressed_genes.txt",
                              header = FALSE,sep = "\t",stringsAsFactors = FALSE)[,1]

RZ_expred_genes_dat <- RZ_g_expr[RZ_expred_genes,-16:-18]

RZ_gene_cor_p <- apply(RZ_expred_genes_dat,1,function(x,y){
  cor.test(x,y)$p.value
},y=rep(c(0,0.01,0.1,1,10),each=3))

RZ_gene_cor_p <- p.adjust(RZ_gene_cor_p,method = "bonferroni")
RZ_gene_cor_sig <- ifelse(RZ_gene_cor_p>0.05,1,2)


RZ_expred_genes_means <- t(apply(RZ_expred_genes_dat,1,function(x){
  tapply(x,INDEX = rep(c("0%","0.01%","0.1%","1%","10%"),each=3),mean)
}))

RZ_expred_genes_means <- as.data.frame(RZ_expred_genes_means,stringsAsFactors=FALSE)

#the palette function would produce an empty Rplots.pdf file in working directory,
#don't know why. If you found a Rplots.pdf file in your working directory, it's not an error
palette(ggplot2::alpha(c("grey","red"),alpha = c(0.05,0.2)))

pdf(paste0(figure_dir,"/RZ_genes_correlated_with_gDNA_trends.pdf"),
    height = 6.5,width = 6.5)
par(mai=par("mai")+0.2)
matplot(t(RZ_expred_genes_means),type = "l",col = RZ_gene_cor_sig,lty = 1,lwd = 2.5,xaxt="n",
        xlab = "Genomic DNA Contamination",ylab = "Gene Expression",font.lab=2,cex.lab=1.5,cex.axis=1.3,cex.main=1.5)
axis(side = 1,labels = c("0%","0.01%","0.1%","1%","10%"),at = axTicks(side = 1),cex.axis=1.3)
legend("topright",legend = c("Correlated with gDNA","Not Correlated with gDNA"),col = 2:1,lwd = 2.5,bty = "n",
       text.font = 2)
dev.off()





#PA genes correlated_with gDNA
PA_pct_expred_genes <- read.table(file = "../../results/expression/PA_pct_only_high_expressed_genes.txt",
                                  header = FALSE,sep = "\t",stringsAsFactors = FALSE)[,1]

PA_pct_expred_genes_dat <- PA_pct_expr[PA_pct_expred_genes,-16:-18]

PA_gene_cor_p <- apply(PA_pct_expred_genes_dat,1,function(x,y){
  cor.test(x,y)$p.value
},y=rep(c(0,0.01,0.1,1,10),each=3))

PA_gene_cor_p <- p.adjust(PA_gene_cor_p,method = "bonferroni")

PA_gene_cor_sig <- ifelse(PA_gene_cor_p>0.05,1,2)

PA_pct_expred_genes_means <- t(apply(PA_pct_expred_genes_dat,1,function(x){
  tapply(x,INDEX = rep(c("0%","0.01%","0.1%","1%","10%"),each=3),mean)
}))

PA_pct_expred_genes_means <- as.data.frame(PA_pct_expred_genes_means,stringsAsFactors=FALSE)

pdf(paste0(figure_dir,"/PA_genes_correlated_with_gDNA_trends.pdf"),
    height = 6.5,width = 6.5)
par(mai=par("mai")+0.2)
matplot(t(PA_pct_expred_genes_means),type = "l",col = PA_gene_cor_sig,lty = 1,lwd = 2.5,xaxt="n",
        xlab = "Genomic DNA Contamination",ylab = "Gene Expression",font.lab=2,cex.lab=1.5,cex.axis=1.3,cex.main=1.5)
axis(side = 1,labels = c("0%","0.01%","0.1%","1%","10%"),at = axTicks(side = 1),cex.axis=1.3)
legend("topright",legend = c("Correlated with gDNA","Not Correlated with gDNA"),col = 2:1,lwd = 2.5,bty = "n",
       text.font = 2)
dev.off()


#expression density plot
pdf(paste0(figure_dir,"/RZ_PA_genes_correlated_with_gDNA_expression_density.pdf"),
    height = 6.5,width = 6.5)
par(mai=par("mai")+0.2)
plot(density(RZ_expred_genes_means[RZ_gene_cor_sig==2,1],na.rm = TRUE),main = "Genes Correlated with gDNA Contamination",
     col=rgb(206,61,50,maxColorValue = 255),lwd=2.5,font.lab=2,cex.lab=1.5,cex.axis=1.3,cex.main=1.5,xlab = NA)
title(xlab = paste0("Ribo-Zero (N = ",sum(RZ_gene_cor_sig==2,na.rm = TRUE),
                    ")\nPoly (A) selection (N = ",sum(PA_gene_cor_sig==2,na.rm = TRUE),")"),
      font.lab=2,cex.lab=1.5,line = 4)
legend("topright",legend = c("Ribo-Zero"),
       col = c(rgb(206,61,50,maxColorValue = 255)),lwd = 2.5,bty = "n",text.font = 2)
dev.off()



## PA and RZ all DEGs
# PA
Poly_A_pct_0_01<-obtain_DEGs_compare(DEG_RData="../../results/expression/DEGs_at_constant_0.01/DEGs_PA_pct.RData",
                                     filter_gene_file="../../results/expression/PA_pct_only_high_expressed_genes.txt")
#RZ
RZ_0_01 <- obtain_DEGs_compare(DEG_RData = "../../results/expression/DEGs_at_constant_0.01/DEGs_RZ.RData",
                               filter_gene_file = "../../results/expression/RZ_only_high_expressed_genes.txt")

pdf(paste0(figure_dir,"/number_of_DEGs.pdf"),width = 6.5,height = 6.5)
DEG_matrix <- matrix(c(sapply(Poly_A_pct_0_01,length),
                       sapply(RZ_0_01,length)),nrow = 2,byrow = TRUE)
clr <- c("#A04EF6","#E7D044")
coor <- barplot(as.vector(DEG_matrix),col = c(rep(clr,4),clr),
                space = c(rep(c(1,0),4),1,0),width = c(rep(1,8)),ylim = c(0,max(DEG_matrix)+1000),
                main = "DEGs between Treatment and Control",cex.main=1.8,sub = "Treatment",
                cex.sub=1.8,font.sub=2,cex.axis = 1.4,
                border = c(rep(clr,4),clr),yaxt="n")
axis(2,at = quantile(0:floor((max(DEG_matrix)+1000)/1000),probs=seq(0,1,0.25))*1000,cex.axis=1.5)
mtext(text = c("0.01%","0.1%","1%","10%","No\nDNase"),side = 1,line = 1.2,
      at = c(2,5,8,11,14),font = 2,cex = 2,padj = 0.5)
legend("topleft",legend = c("Poly (A) selection","Ribo-zero"),
       fill = clr,bty = "n",cex = 2,border = clr)
dev.off()


## DEGs with expression correlated with gDNA or not
PA_cor_genes <- rownames(PA_pct_expred_genes_means)[PA_gene_cor_p<0.05]
RZ_cor_genes <- rownames(RZ_expred_genes_means)[RZ_gene_cor_p<0.05]

write.table(RZ_cor_genes,file = "../../results/expression/RZ_genes_correlated_with_gDNA.txt",col.names = FALSE,
            row.names = FALSE,sep = "\t",quote = FALSE)
write.table(PA_cor_genes,file = "../../results/expression/PA_genes_correlated_with_gDNA.txt",col.names = FALSE,
            row.names = FALSE,sep = "\t",quote = FALSE)


PA_pct_overlap <- lapply(Poly_A_pct_0_01,function(x,y){
  intersect(x,y)
},y=PA_cor_genes)

RZ_overlap <- lapply(RZ_0_01,function(x,y){
  intersect(x,y)
},y=RZ_cor_genes)

PA_pct_unstable <- lapply(Poly_A_pct_0_01,function(x,y){
  setdiff(x,y)
},y=PA_cor_genes)

RZ_unstable <- lapply(RZ_0_01,function(x,y){
  setdiff(x,y)
},y=RZ_cor_genes)

# DEG number
pdf(paste0(figure_dir,"/number_of_DEGs_cor_and_not_cor_with_gDNA_p_adj.pdf"),
    width = 6.5,height = 6.5,colormodel = "cmyk")
DEG_matrix <- matrix(c(sapply(RZ_overlap,length),sapply(RZ_unstable,length)),nrow = 2,byrow = TRUE)

clr <- c("red","grey")
coor <- barplot(DEG_matrix,col = c(rep(clr,4),clr),
                width = c(rep(1,8)),ylim = c(0,max(apply(DEG_matrix,2,sum))+1000),
                main = "DEGs between Treatment and Control",cex.main=1.8,sub = "Treatment",
                #main = "DEGs of Which Expression\nCorrelated with gDNA",cex.main=1.8,sub = "Treatment",
                cex.sub=1.8,font.sub=2,cex.axis = 1.4,
                border = c(rep(clr,4),clr),xaxt="n",yaxt="n")
axis(2,at = quantile(0:floor((max(DEG_matrix)+1000)/1000),probs=seq(0,1,0.25))*1000,cex.axis=1.5)
mtext(text = c("0.01%","0.1%","1%","10%","No\nDNase"),side = 1,line = 1.2,
      at = coor,font = 2,cex = 2,padj = 0.5)

legend("topleft",legend = c("Correlated","Not Correlated"),
       fill = clr,bty = "n",cex = 1.5,border = clr)
dev.off()

# DEG expression density

RZ_overlap <- RZ_overlap[c(-2)]#Drop list element with 1 gene to do density which would cause failure for density

RZ_overlap_expr <- lapply(names(RZ_overlap),function(x,g_expr,DEG_list){
  treatments <- paste0("FPKM.",rep(strsplit(x,split = "_VS_")[[1]],each=3),"_",1:3)
  DEG_exprs <- g_expr[DEG_list[[x]],treatments]
  mean_DEG_exprs <- t(apply(DEG_exprs,1,function(x){
    c(mean(x[1:3]),mean(x[4:6]))
  }))
  colnames(mean_DEG_exprs) <- c("No_gDNA","With_gDNA")
  mean_DEG_exprs
},g_expr=RZ_g_expr,DEG_list=RZ_overlap)
names(RZ_overlap_expr) <- names(RZ_overlap)

RZ_overlap_expr_density <- lapply(RZ_overlap_expr,function(x)list(No_gDNA=density(x[,1]),With_gDNA=density(x[,2])))
RZ_overlap_expr_density_for_compare <- do.call(c,RZ_overlap_expr_density)


#DEGs correlated with gDNA
pdf(paste0(figure_dir,"/RZ_gDNA_cor_DEG_expression_density_p_adj.pdf"),
    width = 13,height = 6.5)
layout(mat = matrix(1:(length(RZ_overlap_expr)+1),nrow = 1,byrow = TRUE),widths = c(0.2,rep(1,3)))
def_mar <- par("mar")
par(mar=rep(0,4))
plot.new()
par(mar=def_mar)
par(mai=c(0.6732,0.15,0.5412,0.1))
plot(RZ_overlap_expr_density[[1]]$No_gDNA,
     ylim = c(0,compare_two_density_y(RZ_overlap_expr_density_for_compare)),
     xlim = c(compare_two_density_x_2(RZ_overlap_expr_density_for_compare[1:2]),
              compare_two_density_x(RZ_overlap_expr_density_for_compare[1:2])),col="blue",
     xlab = NA,ylab = "Density",lwd=2,main = NA,xpd=NA,cex.lab=2.5,cex.axis=1.5,font.lab=2,font.axis=2)
lines(RZ_overlap_expr_density[[1]]$With_gDNA,col="red",lwd=2)
legend("topright",legend = c("0%","0.01%"),title = "gDNA",bty="n",col = c("blue","red"),cex = 2,text.font = 2,lwd = 2)

plot(RZ_overlap_expr_density[[2]]$No_gDNA,
     ylim = c(0,compare_two_density_y(RZ_overlap_expr_density_for_compare)),
     xlim = c(compare_two_density_x_2(RZ_overlap_expr_density_for_compare[3:4]),
              compare_two_density_x(RZ_overlap_expr_density_for_compare[3:4])),col="blue",
     xlab = NA,ylab = NA,lwd=2,main = NA,xpd=NA,cex.lab=2.5,cex.axis=1.5,font.lab=2,font.axis=2,yaxt="n")
lines(RZ_overlap_expr_density[[2]]$With_gDNA,col="red",lwd=2)
legend("topright",legend = c("0%","1%"),title = "gDNA",bty="n",col = c("blue","red"),cex = 2,text.font = 2,lwd = 2)

plot(RZ_overlap_expr_density[[3]]$No_gDNA,
     ylim = c(0,compare_two_density_y(RZ_overlap_expr_density_for_compare)),
     xlim = c(compare_two_density_x_2(RZ_overlap_expr_density_for_compare[5:6]),
              compare_two_density_x(RZ_overlap_expr_density_for_compare[5:6])),col="blue",
     xlab = NA,ylab = NA,lwd=2,main = NA,xpd=NA,cex.lab=2.5,cex.axis=1.5,font.lab=2,font.axis=2,yaxt="n")
lines(RZ_overlap_expr_density[[3]]$With_gDNA,col="red",lwd=2)
axis(side = 2,labels = FALSE)
legend("topright",legend = c("0%","10%"),title = "gDNA",bty="n",col = c("blue","red"),cex = 2,text.font = 2,lwd = 2)

plot(RZ_overlap_expr_density[[4]]$No_gDNA,
     ylim = c(0,compare_two_density_y(RZ_overlap_expr_density_for_compare)),
     xlim = c(compare_two_density_x_2(RZ_overlap_expr_density_for_compare[7:8]),
              compare_two_density_x(RZ_overlap_expr_density_for_compare[7:8])),col="blue",
     xlab = NA,ylab = NA,lwd=2,main = NA,xpd=NA,cex.lab=2.5,cex.axis=1.5,font.lab=2,font.axis=2,yaxt="n")
lines(RZ_overlap_expr_density[[3]]$With_gDNA,col="red",lwd=2)
axis(side = 2,labels = FALSE)
legend("topright",legend = c("0%","No DNase"),title = "gDNA",bty="n",col = c("blue","red"),cex = 2,text.font = 2,lwd = 2)
mtext(text = "Expression",cex=1.5,font=2,side = 1,line = -2,outer = TRUE,adj = c(0.525))

dev.off()

#DEGs not correlated with gDNA
RZ_unstable_expr <- lapply(names(RZ_unstable),function(x,g_expr,DEG_list){
  treatments <- paste0("FPKM.",rep(strsplit(x,split = "_VS_")[[1]],each=3),"_",1:3)
  DEG_exprs <- g_expr[DEG_list[[x]],treatments]
  mean_DEG_exprs <- t(apply(DEG_exprs,1,function(x){
    c(mean(x[1:3]),mean(x[4:6]))
  }))
  colnames(mean_DEG_exprs) <- c("No_gDNA","With_gDNA")
  mean_DEG_exprs
},g_expr=RZ_g_expr,DEG_list=RZ_unstable)
names(RZ_unstable_expr) <- names(RZ_unstable)


RZ_unstable_expr_density <- lapply(RZ_unstable_expr,function(x)list(No_gDNA=density(x[,1]),With_gDNA=density(x[,2])))
RZ_unstable_expr_density_for_compare <- do.call(c,RZ_unstable_expr_density)


pdf(paste0(figure_dir,"/RZ_unstable_DEG_expression_density_p_adj.pdf"),
    width = 13,height = 6.5)

layout(mat = matrix(1:(length(RZ_unstable_expr)+1),nrow = 1,byrow = TRUE),widths = c(0.2,rep(1,5)))
def_mar <- par("mar")
par(mar=rep(0,4))
plot.new()
par(mar=def_mar)
par(mai=c(0.6732,0.15,0.5412,0.1))

plot(RZ_unstable_expr_density[[1]]$No_gDNA,
     ylim = c(0,compare_two_density_y(RZ_unstable_expr_density_for_compare)),
     xlim = c(compare_two_density_x_2(RZ_unstable_expr_density_for_compare[1:2]),
              compare_two_density_x(RZ_unstable_expr_density_for_compare[1:2])),col="blue",
     xlab = NA,ylab = "Density",lwd=2,main = NA,xpd=NA,cex.lab=2.5,cex.axis=1.5,font.lab=2,font.axis=2)
lines(RZ_unstable_expr_density[[1]]$With_gDNA,col="red",lwd=2)
legend("topright",legend = c("0%","0.01%"),title = "gDNA",bty="n",col = c("blue","red"),cex = 2,text.font = 2,lwd = 2)

plot(RZ_unstable_expr_density[[2]]$No_gDNA,
     ylim = c(0,compare_two_density_y(RZ_unstable_expr_density_for_compare)),
     xlim = c(compare_two_density_x_2(RZ_unstable_expr_density_for_compare[3:4]),
              compare_two_density_x(RZ_unstable_expr_density_for_compare[3:4])),col="blue",
     xlab = NA,ylab = NA,lwd=2,main = NA,xpd=NA,cex.lab=2.5,cex.axis=1.5,font.lab=2,font.axis=2,yaxt="n")
lines(RZ_unstable_expr_density[[2]]$With_gDNA,col="red",lwd=2)
axis(side = 2,labels = FALSE)
legend("topright",legend = c("0%","0.1%"),title = "gDNA",bty="n",col = c("blue","red"),cex = 2,text.font = 2,lwd = 2)

plot(RZ_unstable_expr_density[[3]]$No_gDNA,
     ylim = c(0,compare_two_density_y(RZ_unstable_expr_density_for_compare)),
     xlim = c(compare_two_density_x_2(RZ_unstable_expr_density_for_compare[5:6]),
              compare_two_density_x(RZ_unstable_expr_density_for_compare[5:6])),col="blue",
     xlab = "Expression",ylab = NA,lwd=2,main = NA,xpd=NA,cex.lab=2.5,cex.axis=1.5,font.lab=2,font.axis=2,yaxt="n")
lines(RZ_unstable_expr_density[[3]]$With_gDNA,col="red",lwd=2)
axis(side = 2,labels = FALSE)
legend("topright",legend = c("0%","1%"),title = "gDNA",bty="n",col = c("blue","red"),cex = 2,text.font = 2,lwd = 2)

plot(RZ_unstable_expr_density[[4]]$No_gDNA,
     ylim = c(0,compare_two_density_y(RZ_unstable_expr_density_for_compare)),
     xlim = c(compare_two_density_x_2(RZ_unstable_expr_density_for_compare[7:8]),
              compare_two_density_x(RZ_unstable_expr_density_for_compare[7:8])),col="blue",
     xlab = NA,ylab = NA,lwd=2,main = NA,xpd=NA,cex.lab=2.5,cex.axis=1.5,font.lab=2,font.axis=2,yaxt="n")
lines(RZ_unstable_expr_density[[4]]$With_gDNA,col="red",lwd=2)
axis(side = 2,labels = FALSE)
legend("topright",legend = c("0%","10%"),title = "gDNA",bty="n",col = c("blue","red"),cex = 2,text.font = 2,lwd = 2)

plot(RZ_unstable_expr_density[[5]]$No_gDNA,
     ylim = c(0,compare_two_density_y(RZ_unstable_expr_density_for_compare)),
     xlim = c(compare_two_density_x_2(RZ_unstable_expr_density_for_compare[9:10]),
              compare_two_density_x(RZ_unstable_expr_density_for_compare[9:10])),col="blue",
     xlab = NA,ylab = NA,lwd=2,main = NA,xpd=NA,cex.lab=2.5,cex.axis=1.5,font.lab=2,font.axis=2,yaxt="n")
lines(RZ_unstable_expr_density[[5]]$With_gDNA,col="red",lwd=2)
axis(side = 2,labels = FALSE)
legend("topright",legend = c("0%","No DNase"),title = "gDNA",bty="n",col = c("blue","red"),cex = 2,text.font = 2,lwd = 2)
dev.off()





## enriched pathways
warning("The varied version of clusterProfilter would change the number of enriched pathways!")

#DEGs correlated with gDNA enriched pathways
RZ_overlap<-lapply(RZ_0_01,function(x,y){#to produce a new RZ_overlap since we dropped the first 2 elements for density
  intersect(x,y)
},y=RZ_cor_genes)

RZ_overlap_enriched <- lapply(RZ_overlap[1:5],function(x){
  g_id <- gsub("\\..*","",x)
  entrez_id <- try(bitr(geneID = g_id,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")$ENTREZID,
                   silent = TRUE)
  if(class(entrez_id)!="try-error")
  enrichKEGG(gene = entrez_id)
})
RZ_overlap_enriched_pathway_numbers <- sapply(RZ_overlap_enriched,function(x){
  if(is.null(x)){
    0
  }else{
    x <- DOSE:::get_enriched(x)
    length(rownames(x@result))
  }
})

#DEGs not correlated with gDNA enriched pathways
RZ_unstable_enriched <- lapply(RZ_unstable,function(x){
  g_id <- gsub("\\..*","",x)
  entrez_id <- bitr(geneID = g_id,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")$ENTREZID
  
  enrichKEGG(gene = entrez_id)
})
RZ_unstable_enriched_pathway_numbers <- sapply(RZ_unstable_enriched,function(x){
  x <- DOSE:::get_enriched(x)
  length(rownames(x@result))
})

#All DEGs enriched pathways
RZ_0_01_kegg <- lapply(RZ_0_01,function(x){
  g_id <- gsub("\\..*","",x)
  entrez_id <- bitr(geneID = g_id,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")$ENTREZID
  
  enrichKEGG(gene = entrez_id)
})
RZ_0_01_kegg_enriched_pathways_num <- sapply(RZ_0_01_kegg,function(x){
  x <- DOSE:::get_enriched(x)
  length(rownames(x@result))
})

RZ_overall_enriched_pathways <- rbind(RZ_overlap_enriched_pathway_numbers,
                                      RZ_unstable_enriched_pathway_numbers,
                                      RZ_0_01_kegg_enriched_pathways_num)

pdf(paste0(figure_dir,"/Number_of_all_types_DEGs_enriched_pathways_p_adj.pdf"),
    width = 6.5,height = 6.5)
coor <- barplot(RZ_overall_enriched_pathways,col = c("red","grey","blue"),ylim=c(0,max(RZ_overall_enriched_pathways+5)),
                main = "Number of DEG Enriched Pathways",cex.main=1.8,sub = "Treatment",
                cex.sub=1.8,font.sub=2,cex.axis = 1.5,
                border = c("red","grey","blue"),xaxt="n",beside = TRUE)
mtext(text = c("0.01%","0.1%","1%","10%","No\nDNase"),side = 1,line = 1.2,
      at = coor[2,],font = 2,cex = 2,padj = 0.5)
legend("topleft",legend = c("Correlated DEGs","Not Correlated DEGs","All DEGs"),
       fill = c("red","grey","blue"),bty = "n",cex = 1.5,border = c("red","grey","blue"))
dev.off()


system("echo \"cor_genes_and_DEG_number_kegg_and_density.R\" finished! ")