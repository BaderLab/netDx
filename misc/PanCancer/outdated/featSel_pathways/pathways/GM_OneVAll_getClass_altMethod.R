#' alternate version of GM_OneVAll_getClass used for altClassMethod
#' this one doesn't require that the output file for all classes have identical
#' patient sets. 
GM_OneVAll_getClass_alt <- function(resSet) {
    type_rank <- NULL
    for (k in 1:length(resSet)){
    x   <- resSet[[k]]$fullmat
	x <- x[,c("ID","GM_score")]
	colnames(x)[2] <- paste(names(resSet)[k],"SCORE",sep="_")

    if (is.null(type_rank)) {
		type_rank <- x # start building ranking table
	}
    else {
		type_rank <- merge(x=type_rank,y=x,by="ID")
    }
 }

###na_sum <- rowSums(is.na(type_rank[,-1]))
###if (any(na_sum>0))
###    cat(sprintf("*** %i rows have an NA prediction (probably query samples that were not not ranked\n",
###            sum(na_sum>0)))
type_rank <- na.omit(type_rank)

# finally, select the class with the highest rank as the subject label.
maxScore    <- sapply(1:nrow(type_rank),function(i){
                      which.max(type_rank[i,-1])})
patClass    <- sub("_SCORE","",names(maxScore))
type_rank   <- cbind(type_rank, PRED_CLASS=patClass)
type_rank$PRED_CLASS <- as.character(type_rank$PRED_CLASS)
type_rank
}
