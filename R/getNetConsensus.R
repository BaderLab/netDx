#' compile net score across a set of predictor results
#'
#' @details used to compare how individual nets score for different
#' predictor configurations
#' @param scorelist (list) key is dataset name, value is a data.frame 
#' containing PATHWAY_NAME and SCORE. This is the output of
#'  compileFeatureScores()
#' @return (data.frame) Rownames are union of all nets in the input list.
#' Columns show net scores for each key of the input list. Where a 
#' net is not found in a given list, it is assigned the value of NA
#' @examples
#' pathways <- paste("PATHWAY_",1:100,sep="")
#' highrisk <- list()
#' for (k in 1:10) { 
#' 	highrisk[[k]] <- data.frame(PATHWAY_NAME=pathways, 
#'		SCORE=runif(length(pathways),min=0,max=10),
#'				stringsAsFactors=FALSE);
#' }
#' names(highrisk) <- sprintf("Split%i",1:length(highrisk))
#' x <- getNetConsensus(highrisk)
#' @export
getNetConsensus <- function(scorelist) {
    out <- scorelist[[1]]
    colnames(out)[2] <- names(scorelist)[1]
    for (k in 2:length(scorelist)) {
        x <- merge(x = out, y = scorelist[[k]], by = "PATHWAY_NAME", 
						all.x = TRUE, all.y = TRUE)
        colnames(x)[k + 1] <- names(scorelist)[k]
        out <- x
    }
    
    out
}

