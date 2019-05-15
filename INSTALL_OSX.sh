#!/bin/bash

# Script to automate dependency install for netDx

# install R packages
echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
Rscript -e "install.packages(c('devtools','curl','bigmemory','foreach','combinat','doParallel','ROCR','pracma','RColorBrewer','reshape2','ggplot2', 'caroline', 'rmarkdown','igraph','glmnet'))"
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}'
Rscript -e 'if (!requireNamespace("GenomicRanges",quietly=TRUE)) { BiocManager::install("GenomicRanges")}'
Rscript -e 'if (!requireNamespace("RCy3",quietly=TRUE)) { BiocManager::install("RCy3")}'
R CMD INSTALL netDx
R CMD INSTALL netDx.examples
