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
#' @param outPfx (char) prefix for output files. Should include directory
#' name and prefix of file name.
#' @param verbose (logical) print messages
#' @return
#' 1) <outPfx>.gmt file - for enrichment map
#' 2) <outPfx>_nodeAttr.txt (file) table with node properties, notably type,
#' pctPass
#' @examples
#' inDir <- sprintf("%s/extdata/KIRC_output", 
#'	path.package("netDx.examples"))
#' outDir <- paste(getwd(),"plots",sep="/")
#' if (!file.exists(outDir)) dir.create(outDir)
#' featScores <- getFeatureScores(inDir,predClasses=c("SURVIVEYES","SURVIVENO"))
#' gp <- names(featScores)[1]
#' pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt",
#'           path.package("netDx.examples"))
#' pathwayList <- readPathways(pathFile)
#' pathwayList <- pathwayList[c(1:5)]
#' netInfoFile <- sprintf("%s/extdata/KIRC_output/inputNets.txt",
#'      path.package("netDx.examples"))
#' netTypes <- read.delim(netInfoFile,sep="\t",h=FALSE,as.is=TRUE)
#' netInfo <- read.delim(netInfoFile,sep="\t",h=FALSE,as.is=TRUE)
#' output_files <- writeEMapInput(featScores[[gp]],pathwayList,netInfo,
#'                   outPfx=sprintf("%s/%s",outDir,gp))
#' @export
writeEMapInput <- function(featScores, namedSets,netInfo,
	outPfx=sprintf("%s/curr",tempdir()),pctPass=0.70,minScore=1,maxScore=10,
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

	# write node attributes
	netAttrFile <- sprintf("%s_nodeAttrs.txt",outPfx)
	write.table(df2,file=netAttrFile,sep="\t",col=TRUE,row=FALSE,quote=FALSE)

	# write gmt
	outFile <- sprintf("%s.gmt",outPfx)
	conn <- file(outFile)
	for (cur in df2$netName) {
			k2 <- simpleCap(cur)
			if (is.null(namedSets[[cur]])) namedSets[[cur]] <- k2

			curr <- sprintf("%s\t%s\t%s\n", k2,k2,
				paste(namedSets[[cur]],collapse="\t"))
			writeLines(curr,conn)
	}
	close(conn)

  return(c(outFile,netAttrFile))
}

