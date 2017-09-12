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
#' @param pctPass (numeric between 0 and 1) fraction of iterations that
#' a net's score must pass consCutoff for, to be included in the consensus
#' map
#' @param consCutoff (integer) nets must pass this cutoff in at least
#' pctPass of the iterations to be written to file
#' @param maxScore (integer) maximum possible score in one round of cross-
#' validation. e.g. for 10-fold cross-validation, maxScore=10.
#' @param trimFromName (char) strings to trim from name with sub()
#' @param verbose (logical) print messages
#' @return
#' 1) <outPfx>.gmt file - for enrichment map
#' 2) <outPfx>_nodeAttr.txt (file) table with node properties, notably type,
#' pctPass
#' @examples
#' inDir <- sprintf("%s/extdata/KIRC_output", path.package("netDx.examples"))
#' outDir <- paste(getwd(),"plots",sep="/")
#' featScores <- getFeatureScores(inDir,predClasses=c("SURVIVEYES","SURVIVENO"))
#' gp <- names(featScores)[1]
#' pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt",
#'           path.package("netDx.examples"))
#' pathwayList <- readPathways(pathFile)
#' pathwayList <- pathwayList[c(1:5)]
#' netInfo <- read.delim(netInfoFile,sep="\t",h=FALSE,as.is=TRUE)
#' output_files <- writeEMapInput(featScores[[gp]],pathwayList,netInfo,
#'                   outPfx=sprintf("%s/%s",outDir,gp),...)
#' @export
writeEMapInput <- function(featScores, namedSets,netInfo,outPfx="curr",
	pctPass=0.70,maxScore=10,trimFromName=c(".profile","_cont"),verbose=FALSE) {

	dt <- format(Sys.Date(),"%y%m%d")
	netNames <- featScores[,1];
	featScores <- featScores[,-1]

	# compute the max score per net for pctPass % of trials
	maxNetS <- matrix(NA, nrow=length(netNames),ncol=1)
	for (sc in 3:maxScore) {
			tmp <- rowSums(featScores >= sc)
			idx <- which(tmp >= floor(pctPass * ncol(featScores)))
			if (verbose) cat(sprintf("\t%i : %i pass\n", sc, length(idx)))
			maxNetS[idx,1] <- sc
	}
	idx <- which(!is.na(maxNetS))
	maxNetS <- maxNetS[idx,,drop=F]
	netNames <- netNames[idx]

	for (tr in trimFromName) netNames <- sub(tr,"",netNames)

	df1 <- data.frame(netName=netNames, maxScore=maxNetS)
	colnames(netInfo) <- c("netType","netName")
	df2 <- merge(x=df1,y=netInfo,by="netName")

	# write node attributes
	netAttrFile <- sprintf("%s_nodeAttrs_%s.txt",outPfx,dt)
	write.table(df2,file=netAttrFile,sep="\t",col=T,row=F,quote=F)

	# write gmt
	outFile <- sprintf("%s_%s.gmt",outPfx,dt)
	if (file.exists(outFile)) unlink(outFile)
	system(sprintf("touch %s",outFile))
	for (cur in df2$netName) {
			k2 <- .simpleCap(cur)
			if (is.null(namedSets[[cur]])) namedSets[[cur]] <- k2

			cat(sprintf("%s\t%s\t%s\n", k2,k2,
				paste(namedSets[[cur]],collapse="\t")),file=outFile,append=TRUE)
	}
  return(c(outFile,netAttrFile))
}

# convert first letter of each word to uppercase.
 .simpleCap <- function(x) {
				 x <- tolower(x)
         s <- strsplit(x, " ")[[1]]
         paste(toupper(substring(s, 1, 1)), substring(s, 2),
               sep = "", collapse = " ")
     }
