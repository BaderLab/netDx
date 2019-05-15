#!/bin/bash

echo "* Check dependencies"

###if java -version |  grep -q "java version" ; then
###  echo "Java installed."
###else {
###  echo "Install Java before proceeding"
###  exit 0;
### }
###fi
###
###if [[ $(python --version 2>&1) =~ 2\.7 ]]
###    then
###        echo "Python 2.7 installed"
###    else {
###        echo "Install Python 2.7 before proceeding"
###	exit 0;
###	}
###fi
###
###if R -version | grep -q "R version" ; then
###    then
###        echo "R installed"
###    else {
###        echo "Install R before proceeding"
###	exit 0;
###	}
###fi
###
###exit 0

# Script to automate dependency install for netDx

# Must have admin (sudo) access to install the Unix packages
sudo apt-get install zlib1g-dev libssl-dev libssh2-1-dev libcurl-devel

# install R packages
echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
Rscript -e "install.packages(c('devtools','curl','bigmemory','foreach','combinat','doParallel','ROCR','pracma','RColorBrewer','reshape2','ggplot2', 'caroline', 'rmarkdown','igraph','glmnet'))"
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}'
Rscript -e 'if (!requireNamespace("GenomicRanges",quietly=TRUE)) { BiocManager::install("GenomicRanges")}'
Rscript -e 'if (!requireNamespace("RCy3",quietly=TRUE)) { BiocManager::install("RCy3")}'
R CMD INSTALL netDx
R CMD INSTALL netDx.examples
