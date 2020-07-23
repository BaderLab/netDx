# Functions for smoothing mutation networks using interaction networks

#' This function applies the random walk with restart propagation algorithm to a
#' matrix of patients profiles
#'
#' @details A network is an undirected graph G defined by a set of nodes
#'   corresponding to genes, and edges connecting nodes with an experimental
#'   evidence of interaction. A priori nodes are genes for which an information
#'   is known. A novel node is a candidate for being associated to the nodes
#'   above based on their information. A node prediction task leads to detect
#'   novel nodes and propagation techniques are largely applied for the purpose.
#'   Network-based propagation algorithms for node prediction transfer the
#'   information from a priori nodes to any other node in a network. Each node
#'   gets an imputation value which assesses how much information got. The
#'   prediction is based on the guilty-by-association principle. A node with a
#'   high imputation value has a high probability to be associated to a priori
#'   nodes. E.g. in a house where room A has one heater, if room B is the second
#'   hottest room it means that B is close to A and that there is a high
#'   probability that they share a door or wall. These algorithms exploit the
#'   global topology of the network. However, when they are applied to detect if
#'   unknown nodes are functionally associated to known ones, they may suffer of
#'   a drawback depending by the context. In biology, two functionally related
#'   fragments interact physically (direct interaction) or interact indirectly
#'   thanks to one or very few mediators. Therefore, exploring too far
#'   similarities between nodes can introduce noise in the prediction. We apply
#'   a random walk with restart propagation algorithm which resolution is set to
#'   0.2 for giving high values only to the close neighbours of the a priori
#'   nodes.
#' @param mat (data.frame) sparse matrix of patient profiles. Rownames
#'   are unique genes. Colnames are unique patients. A cell is a numeric value.
#' @param net (data.frame) adjancy matrix format of a network
#' @param cl (SOCKcluster) cluster object created with makeCluster function 
#' from parallel
#' @param no_cores (numeric) number of cores used to create the cluster object
#' @return (data.frame) continuous matrix of patient profiles in which each gene
#'   has the final propagation score
#' @importFrom netSmooth netSmooth
#' @rawNamespace import(scater, except = plotHeatmap)
#' @import clusterExperiment
#' @examples 
#'   set.seed(8)
#'   numCores <- 8L
#'   library("netDx")
#'   require("MultiAssayExperiment")
#'   outDir <- "/output"
#'   out_plot_perf <- paste(outDir,"plot_performances.png",sep=getFileSep())
#'   out_features_table <- paste(outDir,"OV_features_table.csv",sep=getFileSep())
#'   out_res_rda <- paste(outDir,"OV_results.rda",sep=getFileSep())
#'   genoFile <- paste(path.package("netDx"),"extdata","OV_mutSmooth_geno.txt",
#'   sep=getFileSep())
#'   geno <- read.delim(genoFile,sep="\t",header=TRUE,as.is=TRUE)
#'   phenoFile <- paste(path.package("netDx"),"extdata","OV_mutSmooth_pheno.txt",
#'   sep=getFileSep())
#'   pheno <- read.delim(phenoFile,sep="\t",header=TRUE,as.is=TRUE)
#'   colnames(geno) <- gsub("-",".",colnames(geno))
#'   pheno$ID <- gsub("-",".",pheno$ID)
#'   rownames(pheno) <- pheno$ID
#'   netFile <- paste(path.package("netDx"),"extdata","CancerNets.txt",
#'   sep=getFileSep())
#'   cancerNets <- read.delim(netFile,sep="\t",header=T,as.is=T)
#'   message("* Excluding genes not present in interaction nets")
#'   noData <- setdiff(rownames(geno),rownames(cancerNets))
#'   if (length(noData)>0) {
#'   	message(paste(length(noData), 
#'   		" genes not present in cancer nets; excluding",sep=""))
#'   	geno <- geno[-which(rownames(geno) %in% noData),]
#'   }
#'   
#'   message("* Running label prop")
#'   require(doParallel)
#'   cl <- makeCluster(numCores)
#'   registerDoParallel(cl)
#'   prop_net <- smoothMutations_LabelProp(geno,cancerNets,cl,no_cores=numCores)
#'   stopCluster(cl)
#' @export
smoothMutations_LabelProp <- function(mat,net,cl,no_cores){
	if (class(mat) == "data.frame") mat <- as.matrix(mat)
	if (class(net) == "data.frame") net <- as.matrix(net)
  #Split the matrix into sections, each one will be processed by one core
  inds <- split(seq_len(ncol(mat)), sort(rep_len(seq_len(no_cores), ncol(mat))))
  res.l <- list()
  #Apply parallelized propagation
  res.l <- foreach(k = 1:length(inds),
	.packages=c("netSmooth","scater","clusterExperiment")) %dopar% {
    nS.res=netSmooth(mat[,inds[[k]]], net , alpha=0.2, verbose = 'auto', 
		normalizeAdjMatrix = c("columns")) 
    return(nS.res)
  }
  #Merge the results
  nS.res=do.call(cbind, res.l)
  #Return the final propagated matrix
  return(nS.res)
}

#' Apply discretization to the matrix resulted from the propagation on the
#' sparse patient matrix
#'
#' @details This function is included in the netDx use case which involves
#'   propagating the sparse matrix of patient's profiles to reduce its sparsity.
#'   This function applies discretization on the propagated matrix of patient
#'   profiles. It sets to 1 the genes which got the highest propagation value.
#'   While, the remaining genes are set to 0. This discretization is driven by
#'   the fact that higher is the propagation value and higher is the chance that
#'   the gene is involved in the patient condition and expression/mutation
#'   profile. On the contrary, genes which got either a medium or a low value
#'   are not trustable.
#' @param m_profiles (data.frame) continous matrix of patient profiles resulted 
#' from applying the network-based propagation algorithm on a binary somatic 
#' mutation sparse matrix.
#' @param m_bin (data.frame) binary somatic mutation sparse matrix. Rownames
#'   are unique genes. Colnames are unique patients. A cell contains a zero or a 
#' one.
#' @return (data.frame) binary somatic mutation matrix which sparsity has been 
#' decreased
#' @examples 
#'   set.seed(8)
#'   numCores <- 8L
#'   library("netDx")
#'   require("MultiAssayExperiment")
#'   outDir <- "/output"
#'   out_plot_perf <- paste(outDir,"plot_performances.png",sep=getFileSep())
#'   out_features_table <- paste(outDir,"OV_features_table.csv",sep=getFileSep())
#'   out_res_rda <- paste(outDir,"OV_results.rda",sep=getFileSep())
#'   genoFile <- paste(path.package("netDx"),"extdata","OV_mutSmooth_geno.txt",
#'   sep=getFileSep())
#'   geno <- read.delim(genoFile,sep="\t",header=TRUE,as.is=TRUE)
#'   phenoFile <- paste(path.package("netDx"),"extdata","OV_mutSmooth_pheno.txt",
#'   sep=getFileSep())
#'   pheno <- read.delim(phenoFile,sep="\t",header=TRUE,as.is=TRUE)
#'   colnames(geno) <- gsub("-",".",colnames(geno))
#'   pheno$ID <- gsub("-",".",pheno$ID)
#'   rownames(pheno) <- pheno$ID
#'   netFile <- paste(path.package("netDx"),"extdata","CancerNets.txt",
#'   sep=getFileSep())
#'   cancerNets <- read.delim(netFile,sep="\t",header=T,as.is=T)
#'   message("* Excluding genes not present in interaction nets")
#'   noData <- setdiff(rownames(geno),rownames(cancerNets))
#'   if (length(noData)>0) {
#'   	message(paste(length(noData), 
#'   		" genes not present in cancer nets; excluding",sep=""))
#'   	geno <- geno[-which(rownames(geno) %in% noData),]
#'   }
#'   
#'   message("* Running label prop")
#'   require(doParallel)
#'   cl <- makeCluster(numCores)
#'   registerDoParallel(cl)
#'   prop_net <- smoothMutations_LabelProp(geno,cancerNets,cl,no_cores=numCores)
#'   stopCluster(cl)
#'   genoP <- thresholdSmoothedMutations(prop_net,geno,"OV_CancerNets")
#' @export
thresholdSmoothedMutations <- function(m_prop,m_bin,name_dataset,n_topXmuts=c(10)){
  m_prop=apply(-m_prop,2,rank)
  n_muts=colSums(m_bin)
  
  m_props_l=list()
  for(k_top in 1:length(n_topXmuts)){
    name_prop=paste(name_dataset,"_x",n_topXmuts[k_top],sep="")
    n_new_muts=n_muts*n_topXmuts[k_top]
    for(i_col in 1:length(n_new_muts)){
      m_prop[m_prop[,i_col]<=n_new_muts[i_col],i_col]=1
      m_prop[m_prop[,i_col]>n_new_muts[i_col],i_col]=0
    }
    m_props_l[[name_prop]]=m_prop
  }
  
  if(length(m_props_l)!=1){
    return(m_props_l)
  }
  if(length(m_props_l)==1){
    return(m_prop)
  }
}