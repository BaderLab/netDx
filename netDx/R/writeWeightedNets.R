#' Write an integrated similarity network consisting of selected networks.
#'
#' @param geneFile (char) path to GENES.txt file created during the making
#' of a generic GeneMANIA database
#' @param netInfo (char) path to NETWORKS.txt file created during the making
#' of a generic GeneMANIA database
#' @param netDir (char) path to directory containing interaction networks.
#' Note that these are networks where the node IDs have been recoded by 
#' GeneMANIA (e.g. 1,2,3)
#' @param keepNets (data.frame) networks to retain in the final network
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
	filterEdgeWt=0,writeAggNet="MAX",writeSingleNets=TRUE,verbose=FALSE){

	pid		<- read.delim(geneFile,sep="\t",header=FALSE,as.is=TRUE)[,1:2]
	colnames(pid)[1:2] <- c("GM_ID","ID")
	netid	<- read.delim(netInfo,sep="\t",header=FALSE,as.is=TRUE)
	colnames(netid)[1:2] <- c("NET_ID", "NETWORK")
	nets	<- merge(x=netid,y=keepNets,by="NETWORK")
	nets	<- nets[,c("NETWORK","NET_ID","WEIGHT")]

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
		nf<- sprintf("%s/1.%i.txt", netDir,nets$NET_ID[i])
		if (verbose) print(basename(nf))
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
				if (verbose) cat(sprintf("%i: %i interactions added\n",
							i, nrow(midx)))
				print(table(maxNet))
				print(summary(as.numeric(intColl)))
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

	# write average PSN
	if (writeAggNet!="NONE") {
		cat("Writing aggregate PSN\n")
		if (writeAggNet=="MEAN") tmp <- intColl/numInt # take mean
		else tmp <- intColl # max value is already in

		require(reshape2)
		tmp[lower.tri(tmp,diag=TRUE)] <- NA # symmetric, remove dups
		ints	<- melt(tmp)
		print(dim(ints))
		cat(sprintf("%i pairs have no interactions\n", 
					sum(is.nan(ints$value))+sum(is.na(ints$value))))
		ints	<- na.omit(ints)
		print(dim(ints))

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
		write.table(ints,file=outF,sep="\t",col=T,row=F,quote=F)
		
		return(outF)
	} else {
		return("")
	}
}

