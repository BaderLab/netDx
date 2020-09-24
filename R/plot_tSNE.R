#' Plot tSNE
#' 
#' @details Plots tSNE of patient similarity network using Rtsne
#' @param psn (matrix) Patient similarity network represented as adjacency
#' matrix (symmetric). Row and column names are patient IDs. Note that NA
#' values will be replaced by very small number (effectively zero).
#' @param pheno (data.frame) Patient labels. ID column is patient ID and 
#' STATUS is patient label of interest. tSNE will colour-code nodes by 
#' patient label.
#' @return (Rtsne) output of Rtsne call. Side effect of tSNE plot
#' @import ggplot2
#' @import Rtsne Rtsne
#' @importFrom RColorBrewer brewer.pal
#' @export
plot_tSNE <- function(psn,pheno) {
message("* Making symmetric matrix")
symmForm <- suppressMessages(makeSymmetric(psn))
symmForm[which(is.na(symmForm))] <- .Machine$double.eps
message("* Running tSNE")
x <- Rtsne(symmForm)
dat <- x$Y
samps <- rownames(symmForm)
idx <- match(samps, pheno$ID)
if (all.equal(pheno$ID[idx],samps)!=TRUE) {
	stop("pheno IDs not matching psn rownames")
}
st <- pheno$STATUS[idx]

# to eliminate the "no visible binding for global variable" problem
y <- status <- NULL

message("* Plotting")
colnames(dat) <- c("x","y")
dat <- as.data.frame(dat,stringsAsFactors=TRUE)
dat$status <- as.factor(st)
p <- ggplot(dat,aes(x,y)) + geom_point(aes(colour=status))
p <- p + xlab("") + ylab("") + ggtitle("Integrated PSN - tSNE")
print(p)

return(x)
}
