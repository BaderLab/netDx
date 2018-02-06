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
sparsify2 <- function(W, outFile="tmp.txt",cutoff=0.3,maxInt=50,EDGE_MAX=1000)  {
	if (maxInt > ncol(W)) maxInt <- ncol(W)

   diag(W) <- 0;
	 W[W < cutoff] <- NA
   x <- apply(W,1,sort,decreasing=TRUE,na.last=NA);
	 for (k in 1:length(x)) {
			cur <- x[[k]]
			tokeep <- names(cur)[1:min(length(cur),maxInt)]
			W[k,which(!colnames(W)%in% tokeep)] <- NA
		}
	tmp <- na.omit(melt(W))
	tmp <- tmp[order(tmp[,3],decreasing=TRUE),]
	#maxEdge <- 0.02*ncol(W); 

	maxEdge <- nrow(tmp)
	if (maxEdge>EDGE_MAX) maxEdge <- EDGE_MAX

	tmp <- tmp[1:maxEdge,]
	write.table(tmp,file=outFile,sep="\t",col=F,row=F,quote=F)

### the code below converts the SIF format back to a matrix,potentially
### for debugging. 
###	W2 <- dcast(tmp,Var2~Var1,value.var="value")
###	rownames(W2) <- W2[,1]; W2 <- W2[,-1]
###	W2 <- W2[,colnames(W)]
###	W2 <- W2[colnames(W),]
###	n <- ncol(W);
###	sp <- nrow(tmp)/(n*(n-1))/2
###	cat(sprintf("%i -> %i edges (%i%% sparsity)\n",
###		sum(!is.na(W)), nrow(tmp), round(sp*100)))
###   return(W2);
}

