#' run lm() for class of interest
#' @param m (data.frame or matrix) rows are measures, columns are samples
#' @param groups (factor) groups to condition model 
#' @param topVar (number 10-100) keep values with top-variance
#' @return output of topTable() for slope. Side effect of plotting pvalue
#' histogram
#' @import limma
runLM <- function(m,groups,topVar=100) {
require(limma)

if (class(groups)=="character") groups <- factor(groups)

if (topVar < 100) {
	if (topVar < 10) topVar <- 10
	cat(sprintf("Limiting top %i %%\n", topVar))
	topVar <- round((topVar/100)*nrow(m))
	mu <- rowMeans(m); n <- ncol(m)
	if (any(is.na(m))) cat("NA detected; removing\n")
	var <- rowSums((m-mu)^2,na.rm=TRUE)/(n-1)
	m <- m[order(var,decreasing=TRUE)[1:topVar],]
}

design <- model.matrix(~groups)
y <- eBayes(lmFit(m,design))
res <- topTable(y,coef=2,num=Inf)

hist(res$P.Value,n=100)
return(res)

}
