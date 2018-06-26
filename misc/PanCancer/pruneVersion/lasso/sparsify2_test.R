#' cleaner sparsification routine
#'
#' @details Sparsifies similarity matrix to keep strongest edges.
#' Sets diagonal and edges < cutoff to NA. Keeps strongest maxInt edges
#' per node. Ties are ignored. Keeps a max of EDGE_MAX edges in the network.
#' @param W (matrix) similarity matrix
#' @param outFile (char) path to file to write sparsified network
#' @param cutoff (numeric) edges with weight smaller than this are set to NA
#' @param maxInt (numeric) max num edges per node.
#' @param EDGE_MAX (numeric) max num edges in network
#' @return writes SIF content to text file (node1,node2,edge weight)
#' @import reshape2
#' @export
sparsify2_test <- function(W, outFile="tmp.txt",cutoff=0.3,maxInt=50,EDGE_MAX=1000)  {
    if (maxInt > ncol(W)) maxInt <- ncol(W)


   diag(W) <- 0;
     W[W < cutoff] <- NA
   x <- apply(W,1,sort,decreasing=TRUE)
	if (x)
browser()
     for (k in 1:length(x)) {
            cur <- x[[k]]
			tryCatch({
            	tokeep <- names(cur)[1:min(length(cur),maxInt)]
			},error=function(ex) { browser()
			})
            W[k,which(!colnames(W)%in% tokeep)] <- NA
        }
    tmp <- na.omit(melt(W))
    tmp <- tmp[order(tmp[,3],decreasing=TRUE),]
    #maxEdge <- 0.02*ncol(W);

    maxEdge <- nrow(tmp)
    if (maxEdge>EDGE_MAX) maxEdge <- EDGE_MAX

    tmp <- tmp[1:maxEdge,]
    write.table(tmp,file=outFile,sep="\t",col=F,row=F,quote=F)
}
