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
#' @param mat (data.frame) Sparse matrix of binarized patient profiles, with
#'	rownames being unique patients and columns, unique genes. Entry [i,j] is
#' 	set to 1 if patient j has a mutation in gene i.
#' @param net (data.frame) Interaction network provided as an adjacency
#' matrix (i.e. symmetric)
#' @param numCores (integer) Number of cores for parallel processing
#' @return (data.frame) Continuous matrix of patient profiles in which each gene
#'   has the final propagation score
#' @importFrom netSmooth netSmooth
#' @rawNamespace import(scater, except = plotHeatmap)
#' @import clusterExperiment
#' @import doParallel
#' @examples 
#' suppressWarnings(suppressMessages(require(MultiAssayExperiment)))
#' require(doParallel)
#' 
#' # load mutation and phenotype data
#' genoFile <- system.file("extdata","TGCT_mutSmooth_geno.txt",package="netDx")
#' geno <- read.delim(genoFile,sep="\t",header=TRUE,as.is=TRUE)
#' phenoFile <- system.file("extdata", "TGCT_mutSmooth_pheno.txt",
#'				package="netDx")
#' pheno <- read.delim(phenoFile,sep="\t",header=TRUE,as.is=TRUE)
#' rownames(pheno) <- pheno$ID
#' 
#' # load interaction nets to smooth over
#' require(BiocFileCache)
#' netFileURL <- paste("http://download.baderlab.org/netDx/",
#' 	"supporting_data/CancerNets.txt",sep="")
#' cache <- rappdirs::user_cache_dir(appname = "netDx")
#' bfc <- BiocFileCache::BiocFileCache(cache,ask=FALSE)
#' rid_rec <- bfcquery(bfc, "CancerNets", "rname")
#' rid <- rid_rec$rid
#' if (!length(rid)) {
#' 	rid <- names(bfcadd(bfc, "CancerNets", netFileURL))
#' }
#' if (!isFALSE(bfcneedsupdate(bfc, rid))){
#' 	bfcdownload(bfc, rid,ask=FALSE)
#' }
#' rid <- rid_rec$rid
#' if (!length(rid)) {
#' 	rid <- names(bfcadd(bfc, "hg18_genes", netFileURL))
#' }
#' if (!isFALSE(bfcneedsupdate(bfc, rid))){
#' 	bfcdownload(bfc, rid,ask=FALSE)
#' }
#' netFile <- bfcrpath(bfc,rids=rid)
#' cancerNets <- read.delim(netFile,sep="\t",header=T,as.is=T)
#' # smooth mutations
#' prop_net <- smoothMutations_LabelProp(geno,cancerNets,numCores=1L)
#' @export
smoothMutations_LabelProp <- function(mat,net,numCores=1L) {
	if (class(mat) == "data.frame") mat <- as.matrix(mat)
	if (class(net) == "data.frame") net <- as.matrix(net)
  #Split the matrix into sections, each one will be processed by one core
  inds <- split(seq_len(ncol(mat)), 
		sort(rep_len(seq_len(numCores), 
		ncol(mat))))

  res.l <- list()

  #Apply parallelized propagation
	cl <- makeCluster(numCores)
	registerDoParallel(cl)

  res.l <- foreach(k = 1:length(inds),
	.packages=c("netSmooth","scater","clusterExperiment")) %dopar% {
    nS.res=netSmooth(mat[,inds[[k]]], net , alpha=0.2, verbose = 'auto', 
		normalizeAdjMatrix = c("columns")) 
    return(nS.res)
  }
	stopCluster(cl)

  #Merge the results
  nS.res <- do.call(cbind, res.l)

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
#' @param smoothedMutProfile (data.frame) continous matrix of patient profiles 
#' resulting from applying :.,$ s/network-based propagation algorithm 
#' (smoothMutations_LabelProp()) on a binary somatic mutation sparse matrix.
#' @param unsmoothedMutProfile (data.frame) binary somatic mutation sparse 
#' matrix. Rownames are unique genes. Colnames are unique patients. A cell 
#' contains a zero or a one.
#' @param nameDataset (char) for titles on plot
#' @param n_topXmuts (numeric between 0 and 1) percent of top mutations
#' to keep. This function converts these to 1.0 when binarizing, so they
#' remain in the thresholded output matrix; other mutations are set to zero.
#' @return (data.frame) binary somatic mutation matrix which sparsity has been 
#' decreased
#' @examples 
#' suppressWarnings(suppressMessages(require(MultiAssayExperiment)))
#' require(doParallel)
#' 
#' # load mutation and phenotype data
#' genoFile <- system.file("extdata","TGCT_mutSmooth_geno.txt",package="netDx")
#' geno <- read.delim(genoFile,sep="\t",header=TRUE,as.is=TRUE)
#' phenoFile <- system.file("extdata", "TGCT_mutSmooth_pheno.txt",
#'				package="netDx")
#' pheno <- read.delim(phenoFile,sep="\t",header=TRUE,as.is=TRUE)
#' rownames(pheno) <- pheno$ID
#' 
#' # load interaction nets to smooth over
#' require(BiocFileCache)
#' netFileURL <- paste("http://download.baderlab.org/netDx/",
#' 	"supporting_data/CancerNets.txt",sep="")
#' cache <- rappdirs::user_cache_dir(appname = "netDx")
#' bfc <- BiocFileCache::BiocFileCache(cache,ask=FALSE)
#' rid_rec <- bfcquery(bfc, "CancerNets", "rname")
#' rid <- rid_rec$rid
#' if (!length(rid)) {
#' 	rid <- names(bfcadd(bfc, "CancerNets", netFileURL))
#' }
#' if (!isFALSE(bfcneedsupdate(bfc, rid))){
#' 	bfcdownload(bfc, rid,ask=FALSE)
#' }
#' rid <- rid_rec$rid
#' if (!length(rid)) {
#' 	rid <- names(bfcadd(bfc, "hg18_genes", netFileURL))
#' }
#' if (!isFALSE(bfcneedsupdate(bfc, rid))){
#' 	bfcdownload(bfc, rid,ask=FALSE)
#' }
#' netFile <- bfcrpath(bfc,rids=rid)
#' cancerNets <- read.delim(netFile,sep="\t",header=T,as.is=T)
#' # smooth mutations
#' prop_net <- smoothMutations_LabelProp(geno,cancerNets,numCores=1L)
#' genoP <- thresholdSmoothedMutations(
#'    prop_net,geno,"TGCT_CancerNets",c(20)
#'   )
#' @export 
thresholdSmoothedMutations <- function(smoothedMutProfile,
		unsmoothedMutProfile,
		nameDataset,n_topXmuts=c(10)){
  smoothedMutProfile=apply(-smoothedMutProfile,2,rank)
  n_muts=colSums(unsmoothedMutProfile)
  
  smoothedMutProfiles_l=list()
  for(k_top in 1:length(n_topXmuts)){
    name_prop=paste(nameDataset,"_x",n_topXmuts[k_top],sep="")
    n_new_muts=n_muts*n_topXmuts[k_top]
    for(i_col in 1:length(n_new_muts)){
      smoothedMutProfile[smoothedMutProfile[,i_col]<=n_new_muts[i_col],i_col]=1
      smoothedMutProfile[smoothedMutProfile[,i_col]>n_new_muts[i_col],i_col]=0
    }
    smoothedMutProfiles_l[[name_prop]]=smoothedMutProfile
  }
  
  if(length(smoothedMutProfiles_l)!=1){
    return(smoothedMutProfiles_l)
  }
  if(length(smoothedMutProfiles_l)==1){
    return(smoothedMutProfile)
  }
}
