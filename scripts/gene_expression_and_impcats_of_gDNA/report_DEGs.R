former_variable <- ls()

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



colnames(all_samples) <- gsub("FPKM\\.","",colnames(all_samples))

##all_samples
all_reps <- substr(colnames(all_samples),1,nchar(colnames(all_samples))-2)

all_sample_names <- colnames(all_samples)

all_sample_DEGs <- dec_all_DEG(expr_dt = all_samples,reps = all_reps,sample_names = all_sample_names)

if(PA_and_RZ){
  save(all_sample_DEGs,file = paste0(figure_folder,"DEGs_",file_name,".RData"))
}else{
  all_DEGS <- all_sample_DEGs[1:5]
  all_DEGS <- lapply(all_DEGS,rownames)
  save(all_DEGS,file = paste0(figure_folder,"DEGs_",file_name,".RData"))
}


all_variable <- ls()
rm(list = c("all_variable",all_variable[!all_variable%in%former_variable]))

