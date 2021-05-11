#' Plots a set of ROC/PR curves with average.
#'
#' @details Plots average curves with individual curves imposed.
#' @param inList (list or ROCR::performance object) ROCR::performance objects, one per iteration
#' @param plotType (char) one of ROC | PR | custom. Affects x/y labels
#' @param xlab (char) x-axis label
#' @param ylab (char) y-axis label
#' @param xlim (numeric) min/max extent for x-axis
#' @param ylim (numeric) min/max extent for y-axis
#' @param plotTitle (numeric) plot title
#' @param meanCol (char) colour for mean trendline
#' @return No value. Side effect of plotting ROC and PR curves
#' @importFrom graphics abline axis par points segments text title hist
#' @examples
#' inDir <- system.file("extdata","example_output",package="netDx")
#' all_rng <- list.files(path = inDir, pattern = 'rng.')
#' fList <- paste(inDir,all_rng,'predictionResults.txt',sep=getFileSep())
#' rocList <- list()
#' for (k in seq_len(length(fList))) {
#'   dat <- read.delim(fList[1],sep='\t',header=TRUE,as.is=TRUE)
#'   predClasses <- c('LumA', 'notLumA')
#'   pred_col1 <- sprintf('%s_SCORE',predClasses[1])
#'   pred_col2 <- sprintf('%s_SCORE',predClasses[2])
#'   idx1 <- which(colnames(dat) == pred_col1)
#'   idx2 <- which(colnames(dat) == pred_col2)
#'  pred <- ROCR::prediction(dat[,idx1]-dat[,idx2], 
#' 		dat$STATUS==predClasses[1])
#'  rocList[[k]] <- ROCR::performance(pred,'tpr','fpr')
#' }
#' plotPerf_multi(rocList,'ROC')
#' @importFrom stats aggregate
#' @importFrom methods slotNames
#' @export
plotPerf_multi <- function(inList, plotTitle = "performance", 
		plotType = "ROC", xlab = "TPR", 
    ylab = "FPR", meanCol = "darkblue", xlim = c(0, 1), ylim = c(0, 1)) {

    if (class(inList)=="performance"){  #if (isa(inList,"performance")){
        inList <- list(inList)
    }
    
    if (plotType == "ROC") {
        xlab <- "TPR"
        ylab <- "FPR"
    } else if (plotType == "PR") {
        xlab <- "Precision"
        ylab <- "Recall"
    } else {
        message("custom type plot\n")
    }
    
    plot(0, 0, type = "n", bty = "n", las = 1, xlim = xlim, ylim = ylim, 
				xlab = xlab, 
        ylab = ylab, main = plotTitle, cex.axis = 1.8,
				cex.lab=1.6)
    out <- list()
    
    is_empty <- 0
    for (k in seq_len(length(inList))) {
        if (length(slotNames(inList[[k]])) == 0) {
            is_empty <- is_empty + 1
            next
        }
        
        x <- inList[[k]]@x.values[[1]]
        y <- inList[[k]]@y.values[[1]]
        
        cur <- aggregate(y, by = list(xvals = x), FUN = mean, na.rm = TRUE)
        colnames(cur) <- c("x", "y")
        out[[k]] <- cur
        
        points(x, y, type = "l", col = "gray90", lwd = 3)
    }
    
    # plot average trendline
    x <- inList[[k]]@x.values[[1]]
    y <- inList[[k]]@y.values[[1]]
    cur_y <- aggregate(y, by = list(xvals = x), FUN = mean, na.rm = TRUE)
    cur <- cbind(cur_y, k)
    colnames(cur) <- c("x", "y", "k")
    out[[k]] <- cur
    points(x, y, type = "l", col = meanCol, lwd = 4)
    
    
    if (length(inList)>1){
    text(0.8 * xlim[2], 0.1 * ylim[2], 
			sprintf("%i splits", length(inList) - is_empty), 
        cex = 1.3)
    }
    
    if (plotType == "ROC") 
        abline(0, 1, col = "red", lwd = 3) else if (plotType == "PR") 
        abline(h = 0.5, col = "red", lwd = 3)
}
