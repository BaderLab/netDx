#' univariate filter by best separation point
#'
#' runs univariate filtering. loops over various fdr cutoffs, evaluating
#' group separation by silhouette score. picks cutoff with best silhouette
#' score
#' @param groups (factor, integer or numeric) clusters
#' @param m (data.frame or matrix) rows are measures, columns are samples
#' @param topVar (number 10-100) keep values with top-variance
#' @return (list) 1) sil_width=avg silh width for various cutoffs
#' 2) bestThresh= cutoff with best silhouette score
#' 3) res: output of runLM().univariate test results 
LMprune <- function(groups,m,topVar=100,verbose=FALSE) {

source("runLM.R")
source("silh.R")
res <- runLM(m,groups,topVar=topVar)

if (min(res$adj.P.Val) > 0.9) {
	cat("smallest FDR is 0.9; not pruning\n")
	return(NA)
}
thresh_vec <- seq(min(res$adj.P.Val),0.9,0.05)
sil_width <- matrix(NA,nrow=length(thresh_vec),ncol=3)
colnames(sil_width) <- c("thresh","num_vars","avg_sil_width")

sil_width[,1] <- thresh_vec
ctr <- 1

xbef <- silh(groups,m,plotMe=FALSE)
ct <- nrow(m)
# evaluate effect of different Q cutoffs
for (thresh in thresh_vec) {
	if (verbose) cat(sprintf("cutoff %1.2f\n", thresh))
	if (sum(res$adj.P.Val < thresh) < 5) {
		sil_width[ctr,2] <- 0
		if (verbose) cat("\t < 5 values left - ignore\n")
	} else {
	res_cur <- subset(res, adj.P.Val < thresh)
	if (verbose) cat(sprintf("\t%i of %i measures left\n",nrow(res_cur), ct,thresh))
	m_cur <- m[which(rownames(m) %in% rownames(res_cur)),]
	x <- silh(groups, m_cur,plotMe=FALSE)
	y <- summary(x)
	if (verbose)cat(sprintf("\tsilh = %1.2f,  %1.2f; avg = %1.2f\n",
			thresh, y$clus.avg.widths[1],y$clus.avg.widths[2],
			y$avg.width))
	sil_width[ctr,2:3] <- c(nrow(m_cur),y$avg.width)
	}
	if (verbose) cat("\n----------------\n")
	ctr <- ctr+1
}

ybef <- summary(xbef)
sil_width <- rbind(c(1,nrow(m),ybef$avg.width),sil_width)

# pick best cutoff
sil_width <- na.omit(sil_width)
idx <- which.max(sil_width[,3])
cat(sprintf("Best cutoff = %1.2f; %i measures; sil width = %1.2f\n",
		sil_width[idx,1], sil_width[idx,2],sil_width[idx,3]))
bestThresh <- sil_width[idx,1]

return(list(sil_width=sil_width,bestThresh=bestThresh,res=res))

}
