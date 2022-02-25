all_args <- commandArgs(trailingOnly = TRUE)

no_merge_gene_fpkm_RData <- all_args[[1]]
figure_dir <- all_args[[2]]

##draw heatmap for RZ and Poly (A) selection
draw_heatmap <- function(expr_dt,reps,batch,n_rep=3,Col_V=FALSE,figure_folder,Row_V=TRUE){
  ##sample orders must 
  ##corresponding to 
  ##the colnames(expr_dt)
  ##make colors
  library(RColorBrewer)

  rep_col <- c(RColorBrewer::brewer.pal(9,"BuPu")[-c(1:4)],"#0868AC")
  
  names(rep_col) <- c("0%","0.01%","0.1%","1%","10%","No DNase")
  ann_col_dat <- data.frame(gDNA=rep(rep(c("0%","0.01%","0.1%","1%","10%","No DNase"),each=3),2),
                        Library=rep(c("Poly (A) selection","Ribo-Zero"),each=18),stringsAsFactors = FALSE)
  rownames(ann_col_dat) <- colnames(expr_dt)
  
  ann_col_col <- list(gDNA=rep_col,Library=setNames(c("#A04EF6","#E7D044"),c("Poly (A) selection","Ribo-Zero")))

  
  ##exclude the row which its sd equal to 0
  kick_out <- sapply(1:nrow(expr_dt),function(x,y){
    if(sd(y[x,])!=0){TRUE}else{message(paste0("Kick the row ",x," out","\n"));FALSE}
  },y=expr_dt)
  
  m0sd1 <- t(apply(expr_dt[kick_out,],1,function(x){
    (x-mean(x))/sd(x)
  }))
  #print(dim(m0sd1))
  library(gplots)
  library(pheatmap)

  pdf(paste0(figure_folder,"heatmap",".pdf"),width = 6,height = 6)

  par(cex.main=2)
  pheatmap(mat = m0sd1,labels_row = rep("",nrow(m0sd1)),labels_col = rep("",ncol(m0sd1)),annotation_col = ann_col_dat[,2:1],
           annotation_row = NA,angle_col = 45,color = colorRampPalette(c("green", "black", "firebrick3"))(101),
           clustering_method = "ward.D",annotation_colors = ann_col_col,breaks=seq(-2,2,length.out=101),
           cluster_cols = Col_V,cluster_rows=Row_V,main = "",border_color = NA,
           cex = 1)
  
  dev.off()
}



#draw PCA with batch
draw_batch_PCA_2 <- function(expr_dt,reps,n_rep=3,bats,figure_folder){
  opar <- par(no.readonly = TRUE)
  #make colors
  # bat_col <- ifelse(grepl(bats[1],colnames(expr_dt)),mk_col(length(unique(bats)),1)[1],
  #                   mk_col(length(unique(bats)),1)[2])
  points_col <- as.integer(factor(reps,levels=unique(reps)))
  palette(c(RColorBrewer::brewer.pal(9,"BuPu")[-c(1:4)],"#0868AC"))
  
  ##calculate pca
  pc.cr<-prcomp(t(expr_dt),retx = TRUE)
  pc<-round(summary(pc.cr)$importance[2,],3)
  
  ##set group
  pt_chr <- sapply(reps,function(x){grep(x,unique(reps))})
  
  
  pdf(paste0(figure_folder,"/PCA",".pdf"),width = 6.5,height = 6.5)
  par(mai=par("mai")+0.2)
  plot(pc.cr$x[,1:2],col=points_col,xlab=paste0("PC1 (",pc[1]*100,"%)"),
       ylab=paste0("PC2 (",pc[2]*100,"%)"),
       pch=rep(c(18,16),each=18),cex.lab=2,font.lab=2,cex=2,cex.axis=1.4,font.axis=2)
  pc_split <- split(pc.cr$x[,1:2],f = reps)
  pc_split <- pc_split[unique(reps)]
  mapply(function(x,y)polygon(x = x[1:n_rep],y=x[((n_rep)+1):(2*n_rep)],border=y,
                              lwd = 1.2),pc_split,unique(points_col))
  legend("topleft",legend = bats,pch = c(18,16),bty="n",cex=2)
  dev.off()
  par(opar)
}

select_genes <- function(expr_dt,expr_pct=0.7,expr_thres=0){
  n_clmn <- ncol(expr_dt)
  #if the expression of one gene in a samples is bigger than expr_thres
  total_count <- apply(expr_dt,1,function(x){
    sum(x > expr_thres)
  })
  total_count[rownames(expr_dt)] > expr_pct*n_clmn
}

dir.create(path = figure_dir,showWarnings = FALSE,recursive = TRUE)

library(ballgown)
load(no_merge_gene_fpkm_RData)
no_merge_gene_fpkm <- log2(no_merge_gene_fpkm+0.01)
colnames(no_merge_gene_fpkm) <- gsub("FPKM.","",colnames(no_merge_gene_fpkm))

all_samples <- no_merge_gene_fpkm

PA_pct_and_RZ_samples <- all_samples[,c(c(28,29,30,1:6,22:24,19:21,25:27),c(46:48,31:36,40:42,37:39,43:45))]

expr_filtered <- PA_pct_and_RZ_samples[select_genes(expr_dt = PA_pct_and_RZ_samples,expr_pct = 0.3,
                              expr_thres = log2(0.02)),]
expr_ordered <- expr_filtered

##HCA
draw_heatmap(expr_dt = expr_ordered,reps = gsub("_[1-3]$","",colnames(expr_ordered)),
             batch = c("Poly (A) selection","Ribo-Zero"),n_rep = 3,Col_V = TRUE,Row_V = TRUE,
             figure_folder = figure_dir)



##PCA label
all_reps <- substr(colnames(expr_ordered),1,nchar(colnames(expr_ordered))-2)
all_bats <- c("Poly (A) selection","Ribo-Zero")


draw_batch_PCA_2(expr_dt = expr_ordered,reps = all_reps,n_rep = 3,bats = all_bats,figure_folder = figure_dir)

system("echo \"PCA_and_heatmap.R finished! \"")
