# convert first letter of each word to uppercase.
 .simpleCap <- function(x) {
				 x <- tolower(x) 
         s <- strsplit(x, " ")[[1]]
         paste(toupper(substring(s, 1, 1)), substring(s, 2),
               sep = "", collapse = " ")
     }	

#' write enrichment map for consensus nets
#' @param netScores (char) path to file with net score)
#' @param setList (list) list of nets and units (e.g.e pathway names and genes)
#' should only contain units profiled in this dataset
#' @param unitList (char) unit list
#' trimFromName (char) strings to trim from name with sub()
#' @return
#' 1) <outPfx>.gmt file - for enrichment map
#' 2) <outPfx>_nodeAttr.txt (file) table with node properties, notably type,
#' pctPass
writeEMap <- function(netScores, setList, netInfo,pctPass=0.70,
		maxScore=10,outPfx="curr",trimFromName=c(".profile","_cont")) {

	dt <- format(Sys.Date(),"%y%m%d")

	netS <- read.delim(netScores,sep="\t",h=T,as.is=T)
	netNames <- netS[,1]; netS <- netS[,-1]
	
	# compute the max score per net for pctPass % of trials
	maxNetS <- matrix(NA, nrow=length(netNames),ncol=1)
	for (sc in 3:maxScore) {
			tmp <- rowSums(netS >= sc)
			idx <- which(tmp >= floor(pctPass * ncol(netS)))
			cat(sprintf("\t%i : %i pass\n", sc, length(idx)))
			maxNetS[idx,1] <- sc
	}
	idx <- which(!is.na(maxNetS))
	maxNetS <- maxNetS[idx,,drop=F]
	netNames <- netNames[idx]

	for (tr in trimFromName) netNames <- sub(tr,"",netNames)

	df1 <- data.frame(netName=netNames, maxScore=maxNetS)
	colnames(netInfo) <- c("netType","netName")
	df2 <- merge(x=df1,y=netInfo,by="netName")
		
	outFile <- sprintf("%s_%s.gmt",outPfx,dt)
	netAttrFile <- sprintf("%s_nodeAttrs_%s.txt",outPfx,dt)

	if (file.exists(outFile)) unlink(outFile)
	system(sprintf("touch %s",outFile))

	# write node attributes
	write.table(df2,file=netAttrFile,sep="\t",col=T,row=F,quote=F)
	# write gmt
	for (cur in df2$netName) {
			k2 <- .simpleCap(cur)
			cat(sprintf("%s\t%s\t%s\n", k2,k2,
				paste(xprList[[cur]],collapse="\t")),file=outFile,append=TRUE)
	}
}


