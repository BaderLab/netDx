#' Write an integrated similarity network consisting of selected networks.
#'
#' @param geneFile (char) path to GENES.txt file created during the making
#' of a generic GeneMANIA database
#' @param netInfo (char) path to NETWORKS.txt file created during the making
#' of a generic GeneMANIA database
#' @param netDir (char) path to directory containing interaction networks.
#' Note that these are networks where the node IDs have been recoded by 
#' GeneMANIA (e.g. 1,2,3)
#' @param keepNets (char or data.frame) networks to include in integrated net
#' If data.frame must be in "NETWORK" column,other columns will be
#' ignored. Mainly included as convenience so pathway scores can passed
#' in table format
#' (NETWORK), and a multiplier constant for edges in that network (WEIGHT)
#' @param outDir (char) path to directory where network files should be 
#' written
#' @param filterEdgeWt (numeric) keep edges with raw edge
#' weight strictly greater than this value. Note that "raw" refers to 
#' this filter being applied before the multiplier is applied.
#' @param writeAggNet (char, one of: [NONE|MEAN|MAX]) Aggregate the network 
#' 1) NONE: does not write aggregate network
#' 2) MEAN: average of weighted edges (raw x netDx score)
#' 3) MAX: max of raw edge weight
#' @param limitToTop (integer) limit to top strongest connections. Set to
#' Inf to list all connections
#' @param writeSingleNets (logical) keep/delete individual recoded nets.
#' TRUE results in each written to its own file; FALSE does not.
#' @param verbose (logical) print messages if TRUE
#' @return If an aggregated network is written (writeAggNet=TRUE),
#' returns the filename of the net. Else returns an empty value.
#' Side effect of writing one tab-delimited file per
#' patient network; file name is <outDir>/<network_name>.txt, and the
#' aggregated net, if requested.
#' File format is:
#' 1) source patient (SOURCE)
#' 2) target patient (TARGET)
#' 3) network name (NET_NAME)
#' 4) weight similarity for the network (WT_SIM)
#' @export
writeWeightedNets <- function(geneFile,netInfo,netDir,keepNets,outDir,
	filterEdgeWt=0,writeAggNet="MAX",limitToTop=50L,
	writeSingleNets=FALSE,verbose=FALSE){

	if (class(keepNets)=="character") {
		keepNets <- data.frame(NETWORK=keepNets,WEIGHT=1)
		keepNets[,1] <- as.character(keepNets[,1])
	}

	pid		<- read.delim(geneFile,sep="\t",header=FALSE,as.is=TRUE)[,1:2]
	colnames(pid)[1:2] <- c("GM_ID","ID")
	netid	<- read.delim(netInfo,sep="\t",header=FALSE,as.is=TRUE)
	colnames(netid)[1:2] <- c("NET_ID", "NETWORK")
	nets	<- merge(x=netid,y=keepNets,by="NETWORK")
	nets	<- nets[,c("NETWORK","NET_ID","WEIGHT")]
	nets$NET_ID <- as.character(nets$NET_ID)

	# clean name
	x <- sub(".profile$","",nets$NETWORK)
	x <- sub("_cont.txt","",x)
	nets$NETWORK_NAME <- x

	numPat	<- nrow(pid) 
	if (!writeAggNet %in% "NONE"){
		# buffer for pairwise interactions
		intColl <- matrix(0,nrow=numPat,ncol=numPat) 
		# num interactions
		numInt <- matrix(0,nrow=numPat,ncol=numPat)
		if (writeAggNet %in% "MAX") {
			intColl <- matrix(NA,nrow=numPat,ncol=numPat)
			maxNet <- matrix(NA,nrow=numPat, ncol=numPat)
		}
	}

	for (i in 1:nrow(nets)) {
		nf<- sprintf("%s/1.%s.txt", netDir,nets$NET_ID[i])
		ints <- read.delim(nf,sep="\t",h=F,as.is=T)
		ints <- subset(ints, ints[,3]>=filterEdgeWt)
		if (nrow(ints)>=1) {
			midx <- rbind(as.matrix(ints[,c(1:2)]),
						  as.matrix(ints[,c(2:1)]))
			if (writeAggNet=="MEAN") {
				# count each edge in both directions so that the upper
				# and lower triangle of the matrix are filled.
				intColl[midx] <- intColl[midx] + ints[,3]#*nets$WEIGHT[i])
				numInt[midx] <- numInt[midx] + 1
			} else if (writeAggNet=="MAX"){
					### cannot run max() like that
				intColl[midx] <- pmax(intColl[midx],ints[,3],na.rm=TRUE)
				numInt[midx] <- numInt[midx] + 1
				maxNet[midx] <- nets$NET_ID[i]
				if (verbose) cat(sprintf("\t%s: %i: %i interactions added\n",
							basename(nf),i, nrow(midx)))
				##print(table(maxNet))
				##print(summary(as.numeric(intColl)))
			}
	
			# resolve to patient name
			midx <- match(ints[,1],pid$GM_ID)
			if (all.equal(pid$GM_ID[midx],ints[,1])!=TRUE) {
				cat("column 1 doesn't match\n")
				browser()
			}
			ints$SOURCE <- pid$ID[midx]; rm(midx)
			
			midx <- match(ints[,2],pid$GM_ID)
			if (all.equal(pid$GM_ID[midx],ints[,2])!=TRUE) {
				cat("column 2 doesn't match\n")
				browser()
			}
			ints$TARGET <- pid$ID[midx]
	
			ints$NETNAME <- nets$NETWORK_NAME[i]
			colnames(ints)[1:2] <- c("NODEID_SOURCE","NODEID_TARGET")
			ints <- ints[,c(4,5,3,6,1,2)]
			colnames(ints)[3:4] <- c("WT_SIM","NETWORK_NAME")
			ints[,3] <- ints[,3]*nets$WEIGHT[i]
	
			if (writeSingleNets) {
			## write output net
			outF <- sprintf("%s/%s_filterEdgeWt%1.2f.txt",outDir,
				nets$NETWORK_NAME[i],filterEdgeWt)
			write.table(ints,file=outF,sep="\t",col=T,row=F,quote=F)
			}
		}
	}
	if (verbose) cat(sprintf("Total of %i nets merged\n", nrow(nets)))

	# write average PSN
	if (writeAggNet!="NONE") {
		cat("\nWriting aggregate PSN\n")
		if (writeAggNet=="MEAN") tmp <- intColl/numInt # take mean
		else tmp <- intColl # max value is already in

		if (limitToTop >= ncol(intColl)) limitToTop <- Inf
		if (!is.infinite(limitToTop)){
			cat(sprintf("* Limiting to top %i edges per patient",
				limitToTop))
			for (k in 1:ncol(intColl)) {
				mytop <- order(tmp[k,],decreasing=TRUE)
				tmp[k,mytop[(limitToTop+1):length(mytop)]] <- NA
			}
		}

		require(reshape2)
		if (is.infinite(limitToTop)) {
			tmp[lower.tri(tmp,diag=TRUE)] <- NA # symmetric, remove dups
		}

		ints	<- melt(tmp)
		cat(sprintf("\n\t%i pairs have no edges\n", 
					sum(is.nan(ints$value))+sum(is.na(ints$value))))
		ints	<- na.omit(ints)

		if (!is.infinite(limitToTop)) {
			x <- paste(ints[,1],ints[,2],sep=".")
			y <- paste(ints[,2],ints[,1],sep=".")
			dup <- intersect(x,y)  # A->B and B->A 
			if (length(dup)>0) {
				cat(sprintf("\tRemoving %i duplicate edges\n",
					length(dup)))
				ints <- ints[-which(y %in% dup),]
			}
		}
		den <- choose(ncol(intColl),2)
		cat(sprintf("\tSparsity = %i/%i (%i %%)\n",
			nrow(ints), den, round((nrow(ints)/den)*100)))

		# resolve to patient name
		midx <- match(ints[,1],pid$GM_ID)
		if (all.equal(pid$GM_ID[midx],ints[,1])!=TRUE) {
			cat("column 1 doesn't match\n")
		}
		ints$SOURCE <- pid$ID[midx]; rm(midx)
		
		midx <- match(ints[,2],pid$GM_ID)
		if (all.equal(pid$GM_ID[midx],ints[,2])!=TRUE) {
			cat("column 2 doesn't match\n")
		}
		ints$TARGET <- pid$ID[midx]
		ints <- ints[,c(4,5,3,1,2)]
		
		outF <- sprintf("%s/aggregateNet_filterEdgeWt%1.2f_%s.txt",
						outDir,filterEdgeWt,writeAggNet)
		if (!is.infinite(limitToTop)) {
			outF <- sub(".txt",sprintf("top%i.txt",limitToTop),outF)
		}
		write.table(ints,file=outF,sep="\t",col=T,row=F,quote=F)
		
		return(outF)
	} else {
		return("")
	}
}

