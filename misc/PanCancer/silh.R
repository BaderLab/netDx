#' silhouette score
#' @param groups (factor, integer or numeric) clusters
#' @param m (data.frame or matrix) rows are measures, columns are samples
#' @param meas (char) dissimilarity measure. [Pearson]
#' @param pal (char) RColorBrewer palette to use
#' @param title (char) silhouette plot title
#' @param plotMe (logical) if true, plots silhouette plot
#' @returns output of cluster::silhouette(). Side effect of plotting
silh <- function(groups,m, meas="Pearson",pal="Dark2",title="",
	plotMe=TRUE) {

if (class(groups)=="numeric") groups <- as.integer(groups)
 if (class(groups)!="factor") groups <- factor(groups)
idx <- order(groups)
groups <- groups[idx]
m <- m[,idx]

if (meas=="Pearson") {
	d <- as.dist(1-cor(na.omit(m)))
} else {
	cat("not implemented\n"); browser()
}

require(RColorBrewer)
si <- silhouette(as.integer(groups), d)
pal <- suppressWarnings(brewer.pal(name=pal,n=length(levels(groups))))
if (plotMe) plot(si,col=pal[groups],main=title)
#print(summary(si))
return(si)

}
