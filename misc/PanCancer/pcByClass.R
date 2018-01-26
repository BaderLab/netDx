#' PCA, top PC projections colour-coded by predicted classes
#' 
#' @param m (matrix or data.frame) data to run PCA on. Rows are measures,
#' columns are samples
#' @param groups (char) value to colour samples by. Should match order of
#' m
#' @return No value. Side effect of 
#' @import RColorBrewer
#' @import combinat
PCbyClass_simple <- function(m,groups,pal="Dark2") {

require(combinat)
require(RColorBrewer)

rs <- rowSums(is.na(m))
cat(sprintf("%i rows have NAs and will be ignored\n", sum(rs>0)))

pr <- prcomp(t(na.omit(m)))
v <- pr$sdev^2
cumvar <- cumsum(v)
cumvar_pct <- cumvar/sum(pr$sdev^2)

par(mfrow=c(2,3))
cb <- combn(4,2)
groups <- factor(groups)
clrs <- suppressWarnings(brewer.pal(name=pal, n=length(levels(groups))))
for (k in 1:ncol(cb)){
	plot(pr$x[,cb[1,k]],pr$x[,cb[2,k]],col=clrs[groups],bty='n',
		 	xlab=sprintf("Proj, PC %i", cb[1,k]),
			ylab=sprintf("Proj, PC %i", cb[2,k]))
	if (k==1) legend("topright", bty='n',fill=clrs,
					legend=levels(groups))
	title(sprintf("PC %i vs PC %i", cb[2,k],cb[1,k]))
}


}
