installed_pkg <- .packages(all.available = TRUE)

bioc_package <- c("ballgown","clusterProfiler")
CRAN_package <- c("RColorBrewer","gplots","pheatmap")

bioc_to_install <- setdiff(bioc_package,installed_pkg)
if(length(bioc_to_install)>0){
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(bioc_to_install)
}

CRAN_to_install <- setdiff(CRAN_package,installed_pkg)
if(length(CRAN_to_install)>0)
  install.packages(c("RColorBrewer","gplots","pheatmap"))
