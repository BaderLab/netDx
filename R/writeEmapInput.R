#' write enrichment map for consensus nets
#'
#' @param featScores (data.frame) network scores across rounds of cross
#' validation. Rows are networks and columns are network name followed by
#' scores for cross-validation rounds. Output of getFeatureScores()
#' @param namedSets (list) list of nets and units (e.g.e pathway names and
#' genes). Should only contain units profiled in this dataset
#' @param netInfo (data.frame) Table of network name (netName) and type
#' (netType). Type is used to assign shapes to nodes:
#'  clinical                                          clinical
#'       rna GUANOSINE_NUCLEOTIDES__I_DE_NOVO__I__BIOSYNTHESIS
#'       rna                              RETINOL_BIOSYNTHESIS
#' @param pctPass (numeric between 0 and 1) fraction of splits for which
#' the highest score for the network is required, for that to be the network's
#' maxScore
#' @param minScore (integer) features with score below this cutoff are
#' excluded from downstream analyses
#' @param maxScore (integer) maximum possible score in one round of cross-
#' validation. e.g. for 10-fold cross-validation, maxScore=10.
#' @param trimFromName (char) strings to trim from name with sub()
#' @param verbose (logical) print messages
#' @examples
#' inDir <- sprintf("%s/extdata/example_output",path.package("netDx"))
#' outDir <- paste(tempdir(),"plots",sep="/")
#' if (!file.exists(outDir)) dir.create(outDir)
#' featScores <- getFeatureScores(inDir,predClasses=c("SURVIVEYES","SURVIVENO"))
#' gp <- names(featScores)[1]
#' pathwayList <- readPathways(getExamplePathways())
#' pathwayList <- pathwayList[c(1:5)]
#' netInfoFile <- sprintf("%s/extdata/example_output/inputNets.txt",
#'      path.package("netDx"))
#' netInfo <- read.delim(netInfoFile,sep="\t",h=FALSE,as.is=TRUE)
#' output_files <- writeEMapInput(featScores[[gp]],pathwayList,netInfo)
#' @export
writeEMapInput <- function(featScores, namedSets,netInfo,pctPass=0.70,
	minScore=1,maxScore=10,
	trimFromName=c(".profile","_cont"),verbose=FALSE) {
	netNames <- featScores[,1];
	featScores <- as.matrix(featScores[,-1])

	# compute the max score per net for pctPass % of trials
	maxNetS <- matrix(NA, nrow=length(netNames),ncol=1)
	for (sc in minScore:maxScore) {
			if (ncol(featScores)>=2) {
				tmp <- rowSums(featScores >= sc)
				idx <- which(tmp >= floor(pctPass * ncol(featScores)))
			} else {
				idx  <- which(featScores >= sc)
			}
			if (verbose) message(sprintf("\t%i : %i pass", sc, length(idx)))
			maxNetS[idx,1] <- sc
	}
	idx <- which(!is.na(maxNetS))
	maxNetS <- maxNetS[idx,,drop=FALSE]
	netNames <- netNames[idx]

	for (tr in trimFromName) netNames <- sub(tr,"",netNames)

	df1 <- data.frame(netName=netNames, maxScore=maxNetS)
	# TODO this colname setting can cause problems when using some pathway lists
	colnames(netInfo) <- c("netType","netName")
	df2 <- merge(x=df1,y=netInfo,by="netName")

	featSet <- list()
	for (cur in df2$netName) {
			k2 <- simpleCap(cur)
			if (is.null(namedSets[[cur]])) namedSets[[cur]] <- k2
			featSet[[k2]] <- namedSets[[cur]]
	}
  out <- list(nodeAttrs=df2,featureSets=featSet)

  return(out)
}

