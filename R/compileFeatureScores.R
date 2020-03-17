#' Tally the score of networks through cross-validation
#'
#' @param fList (char) Vector of paths to GeneMANIA NRANK files
#' @param filter_WtSum (numeric between 5-100) Limit to top-ranked 
#' networks such that cumulative weight is less than this parameter. 
#' e.g. If filter_WtSum=20, first order networks by decreasing weight; 
#' then keep those whose cumulative weight <= 20.
#' @param verbose (logical) print messages
#' @return (data.frame) Feature name and score; includes features that occur
#' at least once in \code{fList}.
#' @examples
#' netDir <- system.file("extdata","GM_NRANK",package="netDx")
#' netFiles <- sprintf('%s/%s', netDir,dir(netDir,pattern='NRANK$'))
#' pTally <- compileFeatureScores(netFiles,verbose=TRUE)
#' print(head(pTally))
#' @export
compileFeatureScores <- function(fList, filter_WtSum = 100, verbose = FALSE) {
    
    if (filter_WtSum < 5) {
        message("filter_WtSum cannot be < 5 ; setting to 5")
        filter_WtSum <- 5
    }
    
    pathwayTally <- list()
    ctr <- 1
    for (fName in fList) {
        tmp <- basename(fName)
        
        try(
					dat <- read.delim(fName, sep = "\t", header = TRUE, 
						as.is = TRUE, skip = 1),silent = TRUE)
        ctr <- ctr + 1
        
        if (!inherits(dat, "try-error")) {
            # file not empty - continue
            if (verbose) {
                message("Net weight distribution:")
                print(summary(dat$Weight))
            }
            
            # actually - it should already be sorted in decreasing 
						# order if we don't reverse
            # it above - but let's sort anyway
            dat <- dat[order(dat$Weight, decreasing = TRUE), ]
            
            cs <- cumsum(dat$Weight)
            keep_max <- which.min(abs(cs - filter_WtSum))
            
            dat <- dat[seq_len(keep_max), ]
            if (verbose) {
                message(sprintf(paste("filter_WtSum = %1.1f; ",
									"%i of %i networks left",sep=""),
									filter_WtSum, nrow(dat), length(cs)))
						}
            
            # put all Network names in pathwaytally. The ones that 
						# are above threshold (Top pathways) get +1
            for (k in dat$Network) {
                if (!k %in% names(pathwayTally)) 
                  pathwayTally[[k]] <- 0
                pathwayTally[[k]] <- pathwayTally[[k]] + 1
            }
            
        }
    }
    out <- unlist(pathwayTally)
    out <- sort(out, decreasing = TRUE)
    out <- data.frame(name = names(out), score = as.integer(out),
		stringsAsFactors=FALSE)
    out[, 2] <- as.integer(as.character(out[, 2]))
    
    out
}
