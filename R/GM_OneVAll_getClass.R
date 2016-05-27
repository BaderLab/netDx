#' assign patient class when ranked by multiple GM predictors
#'
#' @param resSet (list) output of GM_getQueryROC, each key for a different
#' predictor. names(resSet) contain predictor label
#' @return data.frame: ID, GM_score, PRED_CLASS
#' @export
GM_OneVAll_getClass <- function(resSet) {
	type_rank <- NULL
	for (k in 1:length(predRes)){
    x   <- predRes[[k]]$fullmat
    if (is.null(type_rank)) 
        type_rank <- x[,c("ID","GM_score")]
    else {
        if (all.equal(x$ID, type_rank$ID)!=TRUE){ 
            cat("ids don't match"); browser()
        }
        type_rank <- cbind(type_rank, x[,"GM_score"])
    }
        rnkCol <- paste(names(predRes)[k],"SCORE",sep="_") 
        colnames(type_rank)[ncol(type_rank)] <- rnkCol
 }

na_sum <- rowSums(is.na(type_rank[,-1]))
if (any(na_sum>0)) cat(sprintf("*** %i rows have an NA prediction\n",
			sum(na_sum>0)))
type_rank <- na.omit(type_rank)

# finally, select the class with the highest rank as the subject label.
maxScore    <- sapply(1:nrow(type_rank),function(i){
					  which.max(type_rank[i,-1])})
patClass	<- sub("_SCORE","",names(maxScore))
type_rank	<- cbind(type_rank, PRED_CLASS=patClass)
type_rank$PRED_CLASS <- as.character(type_rank$PRED_CLASS)

type_rank
}
