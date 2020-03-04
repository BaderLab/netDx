#' Write an integrated similarity network consisting of selected networks.
#'
#' @param patientID (data.frame) patient identifiers. Columns include
#' internally-generated identifiers (GM_ID) and user-provided identifiers 
#' (ID)
#' @param netInfo (data.frame) network metadata. Columns include
#' internal network name (NET_ID), user-provided name (NETWORK).
#' If a third optional column named "isBinary" is provided, and contains
#' binary values (i.e. 1 and 0), that indicates that the network contains
#' only binary weights and an alternate similarity computation (PropBinary)
#' will be used (see description).
#' @param netDir (char) path to directory containing interaction networks.
#' Note that these are networks where the node IDs have been recoded by 
#' GeneMANIA (e.g. 1,2,3)
#' @param keepNets (char or data.frame) networks to include in integrated net
#' If data.frame must be in "NETWORK" column,other columns will be
#' ignored. Mainly included as convenience so pathway scores can passed
#' in table format
#' (NETWORK), and a multiplier constant for edges in that network (WEIGHT)
#' @param filterEdgeWt (numeric) keep edges with raw edge
#' weight strictly greater than this value. Note that "raw" refers to 
#' this filter being applied before the multiplier is applied.
#' @param aggNetFunc (char, one of: [MEAN|MAX]) Aggregate the network 
#' 2) MEAN: average of weighted edges (raw x netDx score)
#' 3) MAX: max of raw edge weight
#' @param limitToTop (integer) limit to top strongest connections. Set to
#' Inf to list all connections.
#' @param verbose (logical) print messages if TRUE
#' @param plotEdgeDensity (logical) plot density plot of edge weights, one
#' per input net. Used to troubleshoot problems introduced by specific nets.
#' @return (list) 
#' 1) filterEdgeWt (numeric) Value of filterEdgeWt parameter
#' 2) aggNetFunc (char) Value of aggNetFunc parameter
#' 3) limitToTop (integer) Value of limitToTop parameter
#' 4) aggNet (matrix) Value of limitToTop parameter
#' File format is:
#' 1) source patient (SOURCE)
#' 2) target patient (TARGET)
#' 3) network name (NET_NAME)
#' 4) weight similarity for the network (WT_SIM)
#' @importFrom reshape2 melt
#' @export
writeWeightedNets <- function(patientIDs,netIDs,netDir,keepNets,
	filterEdgeWt=0,aggNetFunc="MAX",limitToTop=50L,
	plotEdgeDensity=FALSE,verbose=FALSE){

	pid <- patientIDs
	netid <- netIDs 

	aggNetFunc  <- toupper(aggNetFunc)
 	if (!aggNetFunc %in% c("MEAN","MAX")) {
		stop("aggNetFunc should be one of: MAX|MEAN\n")
	}

	if (missing(keepNets)) {
		keepNets <- netid$NETWORK
	}
	if (class(keepNets)=="character") {
		keepNets <- data.frame(NETWORK=keepNets,WEIGHT=1)
		keepNets[,1] <- as.character(keepNets[,1])
	}

	simMode <- "normal"
	if (any(colnames(netid) %in% "isBinary")) { 
		message("Binary status provided; switching to BinProp mode of similarity!")
		colnames(netid)[3]<- "isBinary"
		simMode <- "BinProp"
	}
	nets	<- merge(x=netid,y=keepNets,by="NETWORK")
	if (simMode =="normal") {
		nets	<- nets[,c("NETWORK","NET_ID","WEIGHT")]
	} else {
		nets	<- nets[,c("NETWORK","NET_ID","WEIGHT","isBinary")]
	}

	nets$NET_ID <- as.character(nets$NET_ID)

	# clean name
	x <- sub(".profile$","",nets$NETWORK)
	x <- sub("_cont.txt","",x)
	nets$NETWORK_NAME <- x

	numPat	<- nrow(pid) 
	# buffer for pairwise interactions
	intColl <- matrix(0,nrow=numPat,ncol=numPat) 
	# num interactions
	numInt <- matrix(0,nrow=numPat,ncol=numPat)

	contNets <- 1:nrow(nets)

	# aggregate nets - binary mode
	if (simMode=="BinProp") {
		intColl <- matrix(0,nrow=numPat,ncol=numPat)
		binNets <- which(nets[,"isBinary"]>0)
		message(sprintf("Got %i binary nets", length(binNets)))
		for (i in binNets) {
			nf<- sprintf("%s/1.%s.txt", netDir,nets$NET_ID[i])
			ints <- read.delim(nf,sep="\t",h=F,as.is=T)
			ints <- subset(ints, ints[,3]>=filterEdgeWt) # probably never needed but
												 		# harmless
			if (nrow(ints)>=1) {
				midx <- rbind(as.matrix(ints[,c(1:2)]),
						 as.matrix(ints[,c(2:1)]))
				intColl[midx] <- intColl[midx] + ints[,3] # increase tally for pairs
														# with shared events
				numInt[midx] <- 1 # count all binary nets exactly once, as we will 
								# create a single measure for them all.
			}
		}

		if (length(binNets)>0) {
			# finally divide by total num binary nets to get proportion similarity
			intColl <- intColl/length(binNets)
			# now convert to the range between filterEdgeWts and 1 to put on
			# par with correlation-based nets
			tmp <- qexp(intColl)
			midx <- which(numInt > 0)
			oldVal <- intColl[midx]
			intColl[midx] <- ((tmp[midx]/max(tmp[midx]))*(1-filterEdgeWt))
			intColl[midx] <- intColl[midx] + filterEdgeWt
		}

		contNets <- setdiff(contNets, which(nets[,"isBinary"]>0))
		message(sprintf("%i continuous nets left",length(contNets)))
	} 

	# aggregate nets - continuous-valued nets
	for (i in contNets) { 
		nf <- sprintf("%s/1.%s.txt", netDir,nets$NET_ID[i])
		ints <- read.delim(nf,sep="\t",h=F,as.is=T)
		oldcount <- nrow(ints)
		ints <- subset(ints, ints[,3]>=filterEdgeWt)
		if (verbose) {
			message(sprintf("Edge wt filter: %i -> %i interactions", 	
				oldcount,nrow(ints)))
		}
		if (nrow(ints)>=1) {
			midx <- rbind(as.matrix(ints[,c(1:2)]),
						  as.matrix(ints[,c(2:1)]))

			if (aggNetFunc=="MEAN") {
				# count each edge in both directions so that the upper
				# and lower triangle of the matrix are filled.
				still_empty <- which(is.na(intColl[midx]))
				if (any(still_empty)) 
					intColl[midx[still_empty]] <- 0
				intColl[midx] <- intColl[midx] + ints[,3]#*nets$WEIGHT[i])
				if (plotEdgeDensity) {
					tmp <- na.omit(as.numeric(ints[,3]))
					hist(tmp,main=nets$NETWORK[i])
				}
				numInt[midx] <- numInt[midx] + 1

			} else if (aggNetFunc=="MAX"){
				intColl[midx] <- pmax(intColl[midx],ints[,3],na.rm=TRUE)
				if (plotEdgeDensity) {
					plot(density(na.omit(as.numeric(ints[,3]))),
						main=nets$NETWORK[i])
				}
				numInt[midx] <- numInt[midx] + 1
				if (verbose) message(sprintf("\t%s: %i: %i interactions added",
							basename(nf),i, nrow(midx)))
			}
		}
	} # finished aggregating nets

	if (verbose) message(sprintf("Total of %i nets merged", nrow(nets)))
	intColl[which(numInt < 1)] <- NA

	# write average PSN
	if (aggNetFunc=="MEAN") tmp <- intColl/numInt # take mean
	else tmp <- intColl # max value is already in

	if (!is.infinite(limitToTop)) {
		if (limitToTop >= ncol(intColl)) limitToTop <- Inf
	}

	# limit top edges
	if (!is.infinite(limitToTop)){
		message(sprintf("* Limiting to top %i edges per patient",
			limitToTop))
		for (k in 1:ncol(intColl)) {
			mytop <- order(tmp[k,],decreasing=TRUE)
			if (limitToTop <= (length(mytop)-1)) {
				tmp[k,mytop[(limitToTop+1):length(mytop)]] <- NA
			}
		}
		}
	
	# now flatten matrix, remove bottom triangle, etc.,
		if (is.infinite(limitToTop)) {
			tmp[lower.tri(tmp,diag=TRUE)] <- NA # symmetric, remove dups
		}
		
		ints	<- melt(tmp)
		message(sprintf("\n\t%i pairs have no edges (counts directed edges)", 
					sum(is.nan(ints$value))+sum(is.na(ints$value))))
		ints	<- na.omit(ints)

		# remove duplicate edges that have been encoded in both directions
		if (!is.infinite(limitToTop)) {
			torm <- c()
			n <- nrow(ints)
			# painfully slow, need way to vectorize this.
			# this pass-through is needed for initial pruning of duplicates
			for (k in 1:(n-1)) {
					dup <- which(ints[(k+1):n,2]==ints[k,1] & 
										ints[(k+1):n,1]==ints[k,2])
					if (any(dup)) torm <- c(torm,dup)	
			}
			
			message(sprintf("\tRemoving %i duplicate edges",
					length(torm)))
			if (length(torm)>0) ints <- ints[-torm,]

		x <- paste(ints[,1],ints[,2],sep=".")
		y <- paste(ints[,2],ints[,1],sep=".")
		z <- which(y %in% x)
		if (any(z)) {
				ints <- ints[-z,]
				message(sprintf("\tSecond pass-through: removed %i more dups", 
					length(z)))
		}	
		x <- paste(ints[,1],ints[,2],sep=".")
		y <- paste(ints[,2],ints[,1],sep=".")
		dup <- intersect(x,y)
			if (length(dup)>0) {
				stop("still have duplicates"); 
			}

		}

		den <- choose(ncol(intColl),2)
		message(sprintf("\tSparsity = %i/%i (%i %%)",
			nrow(ints), den, round((nrow(ints)/den)*100)))

		# resolve to patient name
		midx <- match(ints[,1],pid$GM_ID)
		if (all.equal(pid$GM_ID[midx],ints[,1])!=TRUE) {
			stop("column 1 doesn't match\n")
		}
		ints$SOURCE <- pid$ID[midx]; rm(midx)
		
		midx <- match(ints[,2],pid$GM_ID)
		if (all.equal(pid$GM_ID[midx],ints[,2])!=TRUE) {
			stop("column 2 doesn't match\n")
		}
		ints$TARGET <- pid$ID[midx]
		ints <- ints[,c(4,5,3,1,2)]
		colnames(ints)[1:3] <- c("source","target","weights")
		
		out <- list(
			filterEdgeWt=filterEdgeWt,
			aggNetFunc=aggNetFunc,
			limitToTop=limitToTop,
			aggNet=ints
		)
	return(out)
}

