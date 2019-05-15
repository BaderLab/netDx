# use Dockerized R ("Rocker") as parent image
FROM ubuntu

USER root


ENV DEBIAN_FRONTEND=noninteractive
ENV TZ 'America/New York'

# R pre-requisites
RUN echo $TZ > /etc/timezone && \ 
	apt-get update && \
	apt-get install -y tzdata && \
    apt-get install -y --no-install-recommends \
    gfortran \
	r-base \ 
	openjdk-8-jre \ 
    gcc make g++ \
	zlib1g-dev libssl-dev libssh2-1-dev libcurl4-openssl-dev \
	liblapack-dev liblapack3 libopenblas-base libopenblas-dev \
	libxml2-dev
# && apt-get clean && \
#    rm -rf /var/lib/apt/lists/*


RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages(c('devtools','curl','bigmemory','foreach','combinat','doParallel','ROCR','pracma','RColorBrewer','reshape2','ggplot2', 'caroline', 'rmarkdown'))"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('Biobase','GenomicRanges', 'RCy3'))";
RUN Rscript -e "install.packages(c('pheatmap','RColorBrewer','gProfileR','RJSONIO','httr','ggplot2','glmnet','igraph'))"
RUN Rscript -e "devtools::install_github('cytoscape/r2cytoscape')"
RUN Rscript -e "devtools::install_github('BaderLab/Easycyrest/EasycyRest@0.1')"

# python required for using genemania, pandoc for compiling the html vignettes
RUN apt-get install -y python2.7 python-pip vim git pandoc

# move netDx package and code
WORKDIR /examples
ADD . /examples

# clone the most recent netDx version and install the R package
RUN git clone https://github.com/BaderLab/netDx.git                                                                                                                                                                              
RUN cd netDx && R CMD INSTALL netDx                                                                                                  
RUN cd netDx && R CMD INSTALL netDx.examples