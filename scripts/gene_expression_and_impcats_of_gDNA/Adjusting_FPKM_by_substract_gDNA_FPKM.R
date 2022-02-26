##*****************************
##DEGs from FPKM_ob - FPKM_gDNA
##*****************************

all_args <- commandArgs(trailingOnly = TRUE)

no_merge_gene_fpkm_RData <- all_args[[1]]
figure_dir <- all_args[[2]]


##detect DEG of all sample combinations
dec_all_DEG <- function(expr_dt,reps,sample_names){
  all_cbns <- mkcbns(reps = reps,sample_names = sample_names)
  lst <- list()
  lst_names <- NULL
  for(i in 1:ncol(all_cbns)){
    one_cbn_expr <- expr_of_one_cbn(expr_dt = expr_dt,all_cbns[,i])
    t_test_of_one_cbn <- t_test_for_two_samples(one_cbn_expr[,1:3],one_cbn_expr[,4:6])
    t_test_of_one_cbn$p_adj <- p.adjust(t_test_of_one_cbn$ttestp,method="BH")
    DEG_of_one_cbn <- t_test_of_one_cbn[t_test_of_one_cbn$ttestp<0.05&abs(t_test_of_one_cbn$log2fc)>=1,]
    DEG_of_one_cbn <- DEG_of_one_cbn[order(abs(DEG_of_one_cbn$log2fc),decreasing = TRUE),]
    lst[[i]] <- DEG_of_one_cbn
    lst_names <- c(lst_names,paste0(substr(all_cbns[1,i],1,nchar(all_cbns[1,i])-2),"_VS_",
                                    substr(all_cbns[4,i],1,nchar(all_cbns[4,i])-2)))
    names(lst) <- lst_names
  }
  lst
}

#make combinations
mkcbns <- function(reps,sample_names){
  cbns <- combn(unique(reps),2)
  apply(cbns,2,function(x){
    element1 <- grep(x[1],sample_names,value = TRUE)
    element2 <- grep(x[2],sample_names,value = TRUE)
    c(element1,element2)
  })
}

##obtain expression of one combination
expr_of_one_cbn <- function(expr_dt,sample_names){##the expr_dt contain expression of all samples
  expr_dt[,sample_names]
}

##do t_test for two samples to obtain DEG
t_test_for_two_samples <- function(expr_dt1,expr_dt2){
  bind_dat <- cbind(expr_dt1,expr_dt2)
  #drop genes with sd = 0 across samples
  dropped_rows <- apply(bind_dat,1,sd)==0
  #kept genes that sd is not equal to 0
  bind_dat <- bind_dat[!dropped_rows,]
  expr_dt1 <- expr_dt1[!dropped_rows,]
  expr_dt2 <- expr_dt2[!dropped_rows,]
  
  pvalue <- function(expr1,expr2) {
    var_value <- var.test(expr1,expr2)
    if(sd(expr1)==0&sd(expr2)==0) {var_value$p.value <- 1}
    if(var_value$p.value<0.05){var_equal <- FALSE}else{var_equal <- TRUE}
    obj<-try(t.test(expr1,expr2,var.equal = var_equal,alternative = "two.sided"), silent=TRUE)
    if (is(obj, "try-error")) return(1) else return(obj$p.value)
  }
  ttestfc<-data.frame(genesymbol=rownames(expr_dt1),
                      ttestp=apply(bind_dat,1,
                                   function(x,a,b){
                                     pvalue(x[1:a],x[(a+1):b])
                                   },a=ncol(expr_dt1),b=ncol(bind_dat)),
                      log2fc=rowMeans(expr_dt1)-rowMeans(expr_dt2))
}

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

select_genes <- function(expr_dt,expr_pct=0.7,expr_thres=0){
  n_clmn <- ncol(expr_dt)
  #if the expression of one gene in a samples is bigger than expr_thres
  total_count <- apply(expr_dt,1,function(x){
    sum(x > expr_thres)
  })
  total_count[rownames(expr_dt)] > expr_pct*n_clmn
}

obtain_DEGs_compare_2 <- function(DEG_RData,filter_gene_list){
  load(DEG_RData)
  filter_list <- filter_gene_list
  variables <- ls()
  all_DEGS <- lapply(all_DEGS,function(x){
    intersect(x,filter_list)
  })
  return(sapply(all_DEGS,length))
  rm(list = variables)
}


load(no_merge_gene_fpkm_RData)
RZ_g_expr <- no_merge_gene_fpkm[,31:48]

gDNA_FPKM_value <- read.table("../../results/expression/RZ_gDNA_FPKM_value.txt",header = FALSE,sep = "\t")
RZ_g_expr <- t(apply(RZ_g_expr,1,function(x,y){
  ifelse(x>=y,x-y,0)
},y=gDNA_FPKM_value[,2]))

RZ_g_expr <- log2(RZ_g_expr+0.01)

RZ_g_expr <- RZ_g_expr[,c(16:18,1:6,10:12,7:9,13:15)]
filter_list <- read.table("../../results/expression/RZ_only_high_expressed_genes.txt",
                          header = FALSE,sep = "\t",stringsAsFactors = FALSE)[,1]

#RZ DEGs after DEG adjusting
all_samples <- RZ_g_expr[,1:18]

colnames(all_samples) <- gsub("FPKM\\.","",colnames(all_samples))
##all_samples
all_reps <- substr(colnames(all_samples),1,nchar(colnames(all_samples))-2)
all_sample_names <- colnames(all_samples)
all_sample_DEGs <- dec_all_DEG(expr_dt = all_samples,reps = all_reps,sample_names = all_sample_names)
all_DEGS <- all_sample_DEGs[1:5]
all_DEGS <- lapply(all_DEGS,rownames)
a <- lapply(all_DEGS,function(x,y){
  intersect(x,y)
},y=filter_list)

FPKM_reduced_DEGs <- sapply(a,length)

#RZ unreduced DEGs
RZ_0_01 <- obtain_DEGs_compare(DEG_RData = "../../results/expression/DEGs_at_constant_0.01/DEGs_RZ.RData",
                               filter_gene_file = "../../results/expression/RZ_only_high_expressed_genes.txt")

FPKM_unreduced_DEGs <- sapply(RZ_0_01,length)

pdf(paste0(figure_dir,"/RZ_DEG_numbers_before_and_after_adjusting.pdf"),width = 6.5,height = 6.5)
clr <- RColorBrewer::brewer.pal(9,"Set1")[1:2]
DEG_matrix <- rbind(FPKM_unreduced_DEGs,FPKM_reduced_DEGs)
coor <- barplot(as.vector(DEG_matrix),col = c(rep(clr,4),clr),
                space = c(rep(c(1,0),4),1,0),width = c(rep(1,8)),ylim = c(0,max(DEG_matrix)+1000),
                main = "DEGs between Treatment and Control",cex.main=1.8,sub = "Treatment",
                cex.sub=1.8,font.sub=2,cex.axis = 1.4,
                border = c(rep(clr,4),clr),yaxt="n")
axis(2,at = quantile(0:floor((max(DEG_matrix)+1000)/1000),probs=seq(0,1,0.25))*1000,cex.axis=1.5)
ylim <- c(0,max(DEG_matrix)+1000)
mtext(text = c("0.01%","0.1%","1%","10%","No\nDNase"),side = 1,line = 1.2,
      at = c(2,5,8,11,14),font = 2,cex = 2,padj = 0.5)

legend("topleft",legend = c("FPKM No Adjusting","FPKM Adjusted"),
       fill = clr,bty = "n",cex = 2,border = clr)
dev.off()





PA_g_expr <-  no_merge_gene_fpkm[,1:30][,c(28:30,1:6,22:24,19:21,25:27)]

gDNA_FPKM_value <- read.table("../../results/expression/PA_gDNA_FPKM_value.txt",header = FALSE,sep = "\t")
gDNA_FPKM_value <- gDNA_FPKM_value[c(16:18,1:6,10:12,7:9,13:15),]
PA_g_expr <- t(apply(PA_g_expr,1,function(x,y){
  ifelse(x>=y,x-y,0)
},y=gDNA_FPKM_value[,2]))

PA_g_expr <- log2(PA_g_expr+0.01)

#PA_g_expr <- PA_g_expr[,c(16:18,1:6,10:12,7:9,13:15)]
filter_list <- read.table("../../results/expression/PA_pct_only_high_expressed_genes.txt",
                          header = FALSE,sep = "\t",stringsAsFactors = FALSE)[,1]


#RZ DEGs after DEG adjusting
all_samples <- PA_g_expr[,1:18]

colnames(all_samples) <- gsub("FPKM\\.","",colnames(all_samples))
##all_samples
all_reps <- substr(colnames(all_samples),1,nchar(colnames(all_samples))-2)
all_sample_names <- colnames(all_samples)
all_sample_DEGs <- dec_all_DEG(expr_dt = all_samples,reps = all_reps,sample_names = all_sample_names)
all_DEGS <- all_sample_DEGs[1:5]
all_DEGS <- lapply(all_DEGS,rownames)
a <- lapply(all_DEGS,function(x,y){
  intersect(x,y)
},y=filter_list)

FPKM_reduced_DEGs <- sapply(a,length)


#PA unreduced DEGs
PA_pct_0_01 <- obtain_DEGs_compare(DEG_RData = "../../results/expression/DEGs_at_constant_0.01/DEGs_PA_pct.RData",
                                   filter_gene_file = paste0("../../results/expression/",
                                                             "PA_pct_only_high_expressed_genes.txt"))
FPKM_unreduced_DEGs <- sapply(PA_pct_0_01,length)

pdf(paste0(figure_dir,"/PA_DEG_numbers_before_and_after_adjusting.pdf"),width = 6.5,height = 6.5)
clr <- RColorBrewer::brewer.pal(9,"Set1")[1:2]
DEG_matrix <- rbind(FPKM_unreduced_DEGs,FPKM_reduced_DEGs)
coor <- barplot(as.vector(DEG_matrix),col = c(rep(clr,4),clr),
                space = c(rep(c(1,0),4),1,0),width = c(rep(1,8)),ylim = ylim,
                main = "DEGs between Treatment and Control",cex.main=1.8,sub = "Treatment",
                cex.sub=1.8,font.sub=2,cex.axis = 1.4,
                border = c(rep(clr,4),clr),yaxt="n")
axis(2,at = quantile(0:floor(max(ylim)/1000),probs=seq(0,1,0.25))*1000,cex.axis=1.5)
mtext(text = c("0.01%","0.1%","1%","10%","No\nDNase"),side = 1,line = 1.2,
      at = c(2,5,8,11,14),font = 2,cex = 2,padj = 0.5)

legend("topleft",legend = c("FPKM No Adjusting","FPKM Adjusted"),
       fill = clr,bty = "n",cex = 2,border = clr)
dev.off()



##************************
##reduce DEGs by threshold
##************************

load(no_merge_gene_fpkm_RData)
no_merge_gene_fpkm <- log2(no_merge_gene_fpkm+0.01)

RZ_g_expr <- no_merge_gene_fpkm[,31:48][,c(16:18,1:6,10:12,7:9,13:15)]
PA_g_expr <- no_merge_gene_fpkm[,1:30][,c(28:30,1:6,22:24,19:21,25:27)]

thres <- seq(0,1,0.01)

##Poly (A) Selection
PA_pct_DEGs <- list()
for(threshold in thres){
  considered_gene <- select_genes(expr_dt = PA_g_expr,expr_pct = 0.3,expr_thres = log2(threshold+0.01))
  filter_list <- rownames(PA_g_expr)[considered_gene]
  
  res <- obtain_DEGs_compare_2(DEG_RData = "../../results/expression/DEGs_at_constant_0.01/DEGs_PA_pct.RData",
                               filter_gene_list = filter_list)
  
  PA_pct_DEGs <- c(PA_pct_DEGs,list(res))
}
names(PA_pct_DEGs) <- thres
PA_pct_mat <- sapply(PA_pct_DEGs,function(x)x)

##Ribo-Zero 
RZ_g_DEGs <- list()
for(threshold in thres){
  considered_gene <- select_genes(expr_dt = RZ_g_expr,expr_pct = 0.3,expr_thres = log2(threshold+0.01))
  filter_list <- rownames(RZ_g_expr)[considered_gene]
  
  res <- obtain_DEGs_compare_2(DEG_RData = "../../results/expression/DEGs_at_constant_0.01/DEGs_RZ.RData",
                               filter_gene_list = filter_list)
  
  RZ_g_DEGs <- c(RZ_g_DEGs,list(res))
}
names(RZ_g_DEGs) <- thres
RZ_g_mat <- sapply(RZ_g_DEGs,function(x)x)


#plot
palette(c(RColorBrewer::brewer.pal(9,"BuPu")[-c(1:4)],"#0868AC","#D23C37","#92D050"))
pdf(paste0(figure_dir,"/PA_pct_DEGs_at_different_thres.pdf"),width = 6.5,height = 6.5)
par(mai=par("mai")+0.2)
plot(y = PA_pct_mat[4,],x = thres,col=5,pch=16,type = "l",lwd=2.5,xlab = "Filter Threshold",
     ylab = "Number of DEGs",font = 2,
     font.lab = 2,cex.main=1.5,ylim=c(0,max(PA_pct_mat)),cex.lab=1.5,cex.axis = 1.3,yaxt="n")
mapply(function(x,y)points(x = thres,y = x,col=y,pch=16,type = "l",lwd=2.5),
       as.data.frame(t(PA_pct_mat[-4,]),stringsAsFactors = FALSE),c(2:4,6))
axis(2,at = quantile(0:floor((max(PA_pct_mat))/100),probs=seq(0,1,0.25))*100,cex.axis=1.5,font.axis=2)
legend("topright",legend = c(paste0(c("0.01%","0.1%","1%","10%")," gDNA"),"No DNase Treatment"),
       col = 2:6,lwd = 2.5,bty="n",text.font = 2)
dev.off()

pdf(paste0(figure_dir,"/RZ_DEGs_at_different_thres.pdf"),width = 6.5,height = 6.5)
par(mai=par("mai")+0.2)
plot(y = RZ_g_mat[4,],x = thres,col=5,pch=16,type = "l",lwd=2.5,xlab = "Filter Threshold",
     ylab = "Number of DEGs",font = 2,
     font.lab = 2,cex.main=1.5,ylim=c(0,max(RZ_g_mat)),cex.lab=1.5,cex.axis = 1.3,yaxt="n")
mapply(function(x,y)points(x = thres,y = x,col=y,pch=16,type = "l",lwd=2.5),
       as.data.frame(t(RZ_g_mat[-4,]),stringsAsFactors = FALSE),c(2:4,6))
axis(2,at = quantile(0:floor((max(RZ_g_mat))/1000),probs=seq(0,1,0.25))*1000,cex.axis=1.5,font.axis=2)
legend("topright",legend = c(paste0(c("0.01%","0.1%","1%","10%")," gDNA"),"No DNase Treatment"),
       col = 2:6,lwd = 2.5,bty="n",text.font = 2)
dev.off()





system("echo \"Adjusting_FPKM_by_substract_gDNA_FPKM.R\" finished! ")

