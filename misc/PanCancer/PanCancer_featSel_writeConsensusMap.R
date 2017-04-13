#' write an enrichment map of the nets that on average score above
#' the cutoff, across a set of iterations.
#'
require(netDx)
require(netDx.examples)
#' @param datDir (char) path to output results directory. Directory structure
#' should be:  <datDir>/rng<runNumber>/<className>/GM_results/*_pathway_CV_score.txt
#' @param predClasses (char) classes that are being predicted
#' @param consCutoff (integer) nets must pass this cutoff in at least 
#' pctPass of the iterations to be written to file
#' @param pctPass (numeric between 0 and 1) fraction of iterations that
#' a net's score must pass consCutoff for, to be included in the consensus
#' map
#' @param netInfo (list)
#' keys are user-defined network types (e.g. clinical, mutations, expression),
#' value is a list with the entry:
#' 1) pattern: (char) a unique pattern to grep networks of this type in 
#' pathway_CV_score.txt. See note with ***.
#' 2) netSet: (list or char) key-value pair of network names and included units. 
#' e.g. pathway names and member genes.
#' If set to "base", the network name is written as is to the output file
#' *** Order in this list matters. If a net category does not have a unique
#' pattern, list it after others that do. Then those nets will be first
#' accounted for, and the residuals will be grep-ed for with this pattern.
#' For example, if mutation nets are MUT_*_cont and clinical nets have *_cont
#' then list the mutation ones first. After these are removed, any with *_cont
#' remaining should be clinical.
#' @param outPfx (char) prefix to output filename (absolute path)
#' @param trimFromName (char) strings to trim from net names before writing
#' to file
#' @return no values. Creates two output files for SURVIVE
writeConsensusMap <- function(datDir,predClasses=NULL,consCutoff=7L, 
	pctPass=0.50,netInfo,outPfx="./pred_",trimFromName=c(".profile","_cont")) {

# collect nets for each group
for (gp in predClasses) {
cat("----------------------\n")
cat(sprintf("%s\n",gp))
cat("----------------------\n")
	rngDirs <- dir(path=datDir, pattern="^rng")
	cat(sprintf("Got %i iterations\n", length(rngDirs)))

	# collect nets from all iterations
	cat("* Collecting single runs\n")
	netColl <- list()
	for (curDir in rngDirs) {
		scoreFile <- sprintf("%s/%s/%s/GM_results/%s_pathway_CV_score.txt",
				 datDir,curDir,gp,gp)
		tmp	 <- read.delim(scoreFile,sep="\t",h=T,as.is=T)
		colnames(tmp)[1] <- "PATHWAY_NAME"
		netColl[[curDir]] <- tmp
	}

	# filter for nets meeting cutoff criteria
	cat("* Computing consensus\n")
	cons <- getNetConsensus(netColl); x1 <- nrow(cons)
	na_sum <- rowSums(is.na(cons))
	cons <- cons[which(na_sum < 1),]
full_cons <- cons
	cat(sprintf("\t%i of %i scored in all rounds\n",length(cons),x1))
	
	avg_score <- rowMeans(cons[,-1])
	passes_cutoff <- cons[,-1]>= consCutoff;
	idx <- which(rowSums(passes_cutoff)>=round(pctPass*length(rngDirs)))
	cat(sprintf("\tConsensus: %i nets scored >= %i\n",length(idx),consCutoff))
	cons <- cons[idx,1]

	for (curT in trimFromName) cons <- sub(curT,"",cons)

	outFile <- sprintf("%s_%s_consensusNets_cutoff%i.gmt",
										 outPfx,gp,consCutoff)
	write.table(cons,file=outFile,sep="\t",col=F,row=F,quote=F)
	netTypeFile <- sub(".gmt",".netTypes.txt",outFile)

	if (file.exists(outFile)) unlink(outFile)
	system(sprintf("touch %s",outFile))
	if (file.exists(netTypeFile)) unlink(netTypeFile)
	system(sprintf("touch %s",netTypeFile))

	# convert first letter of each word to uppercase.
  .simpleCap <- function(x) {
				 x <- tolower(x) 
         s <- strsplit(x, " ")[[1]]
         paste(toupper(substring(s, 1, 1)), substring(s, 2),
               sep = "", collapse = " ")
     }	

	# first write RNA-based pathways
	unaccounted <- cons
	for (curType in names(netInfo)) {
		curDat <- netInfo[[curType]]
		idx <- grep(curDat$pattern,unaccounted)
		cat(sprintf("\t%s (pattern=%s) : %i nets ->", 
								curType,curDat$pattern,length(idx)))

		curNets <- unaccounted[idx]
		if (class(curDat$netSet)=="character") {
				cat("writing single\n")
			for (k in curNets) {
				k2 <- .simpleCap(k)
				cat(sprintf("%s\t%s\t%s\n", k2,k2,k2),file=outFile,append=TRUE)
				cat(sprintf("%s\t%s\n",k2,curType),file=netTypeFile,append=TRUE)
			}
		} else {
			midx <- which(names(curDat$netSet) %in% curNets)
			netSet <- curDat$netSet[midx]
			cat(sprintf("%i matches found\n", length(midx)))
			for (k in names(netSet)) {
				k2 <- .simpleCap(k)
				cat(sprintf("%s\t%s\t%s\n", k2,k2,paste(netSet[[k]],collapse="\t")),
					file=outFile,append=TRUE)
				cat(sprintf("%s\t%s\n",k2, curType),file=netTypeFile,append=TRUE)
			}
		}
		unaccounted <- setdiff(unaccounted, unaccounted[idx])
		cat(sprintf("\t\t...%i still unaccounted for\n\n",length(unaccounted)))
	}
}
}

