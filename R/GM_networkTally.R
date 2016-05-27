#' Tally the score of networks through cross-validation
#'
#' @param fList (char) Vector of paths to GeneMANIA NRANK files
#' @param filter_WtSum (numeric between 5-100) Limit to top-ranked 
#' networks such that cumulative weight is less that this parameter. 
#' e.g. If filter_WtSum=20, first order networks by decreasing weight; 
#' then keep those whose cumulative weight <= 20.
#' @param verbose (logical) print messages
#' @return (named integer) Vector of scores for networks that occur at 
#' least once in \code{fList}.
#' @export
GM_networkTally <- function(fList,filter_WtSum=100,verbose=FALSE) {

if (filter_WtSum < 5) {
	cat("filter_WtSum cannot be < 5 ; setting to 5\n")
	filter_WtSum <- 5;
}
	
pathwayTally <- list()
ctr <- 1
for (fName in fList) {
	tmp	<- basename(fName)

	dat <- try(read.delim(fName,sep="\t",h=T,as.is=T),silent=TRUE)
	ctr <- ctr+1

#print(summary(dat))
	
	if (!inherits(dat,"try-error")) { # file not empty
		# remove first group related line
		dat <- dat[-1,]
		dat <- dat[order(dat$Weight,decreasing=TRUE),]
	
		cs			<- cumsum(dat$Weight)
		keep_max	<- which.min(abs(cs-filter_WtSum))
		dat			<- dat[1:keep_max,]
		if (verbose) cat(sprintf(
				"filter_WtSum = %1.1f; %i of %i networks left\n",
				filter_WtSum, nrow(dat),length(cs)))
		
		for (k in dat[,2]) {
			if (!k %in% names(pathwayTally)) pathwayTally[[k]] <- 0;
			pathwayTally[[k]]<- pathwayTally[[k]]+1;
		}
	}
}

out <- unlist(pathwayTally)
out <- sort(out,decreasing=TRUE)
out <- data.frame(name=names(out), score=as.integer(out))
out[,2] <- as.integer(as.character(out[,2]))

out
}
