#' generate a psn of cosine similarity by concatenating consensus nets
#' for the cancer.

rm(list=ls())
# --------------------------------------------------------------



# --------------------------------------------------------------
# Param for computing integrated PSN
outDir		<- "."
cutoff 		<-10  	# include nets with score >= this value
corrThresh 	<-0.7 	# exclude edges with similarity < this threshold
dt <- format(Sys.Date(),"%y%m%d")
verbose <- FALSE

# --------------------------------------------------------------
# setup for network generation in Cytoscape
require(EasycyRest) # https://github.com/shraddhapai/EasycyRest
require(httr)

styleName <- "PSNstyle"
portNum <- 1234
base.url <- sprintf("http://localhost:%i/v1",portNum)
res	<- httr::GET(sprintf("%s/styles",base.url))
curStyles <- gsub("\\\"","",rawToChar(res$content))
curStyles <- unlist(strsplit(curStyles,","))

if (any(grep("PSNstyle",curStyles))) {
	cat("style exists not creating\n")
} else {
	cat("Creating style\n")
	nodeFills <- map_NodeFillDiscrete("GROUP",
		c("SURVIVEYES","SURVIVENO"),
		c("#FF9900","#003399"))
	defaults <- list("NODE_SHAPE"="ellipse",
			"NODE_SIZE"=50,
			"EDGE_TRANSPARENCY"=120,
			"NODE_TRANSPARENCY"=120)
	sty <- createStyle(styleName, 
		defaults=defaults,
		mappings=list(nodeFills))
}

# --------------------------------------------------------------
# data dirs for input
rootDir <- "/Users/shraddhapai/Documents/Research/BaderLab"
# to create pheno table
clinList <- list(
	KIRC=sprintf("%s/2017_TCGA_KIRC/input/KIRC_clinical_core.txt",rootDir),
	OV=sprintf("%s/2017_TCGA_OV/input/OV_clinical_core.txt",rootDir),
	LUSC=sprintf("%s/2017_TCGA_LUSC/input/LUSC_clinical_core.txt",rootDir),
	GBM=sprintf("%s/2017_TCGA_GBM/input/GBM_clinical_core.txt",rootDir)
)
survList <- list(
	KIRC=sprintf("%s/2017_TCGA_KIRC/input/KIRC_binary_survival.txt",rootDir),
	OV=sprintf("%s/2017_TCGA_OV/input/OV_binary_survival.txt",rootDir),
	LUSC=sprintf("%s/2017_TCGA_LUSC/input/LUSC_binary_survival.txt",rootDir),
	GBM=sprintf("%s/2017_TCGA_GBM/input/GBM_binary_survival.txt",rootDir)
)
# best performing iteration.
selIter <- list(
	KIRC=sprintf("%s/2017_TCGA_KIRC/output/featSel_170222/rng29", rootDir),
	OV=sprintf("%s/2017_TCGA_OV/output/featSel_170327/rng57", rootDir),
	LUSC=sprintf("%s/2017_TCGA_LUSC/output/featSel_incMutRPPA_round2170327/rng4",
							 rootDir),
	GBM=sprintf("%s/2017_TCGA_GBM/output/featSel_incMut_round2_170223/rng78", 
							 rootDir)
)

require(netDx)

# --------------------------------------------------------------
# Work begins

logFile <- sprintf("%s/getPSN_%s.log",outDir,dt)
sink(logFile,split=TRUE)

tryCatch({
# aggFun aggregate edges between any given pair of patients
# simMode normal|BinProp. BinProp converges all binary sims into
		# "proportion binary with similarity".

datSets <- c("OV","KIRC","LUSC","GBM")
dijk <- list()
for (aggFun in c("dummy")) {
	dijk[[aggFun]] <- list()
	for (simMode in "dummy") { 		
		curDijk <- matrix(NA,nrow=length(datSets),ncol=8)
		rownames(curDijk) <- datSets
		colnames(curDijk) <- c("NO","YES","YES-NO","overall","pYES-YESNO",
			"pYES-overall","pNO-YESNO","pNO-overall")

		cur_i <- 1

		for (curSet in datSets) {
		clinical_file <- clinList[[curSet]]
		survival_file <- survList[[curSet]]
		dataDir <- selIter[[curSet]]

		cat(sprintf("%s: %s: %s\n", curSet, simMode, aggFun))

		# -----------------------------------------------------
		# patient IDs - should be identical for both
		ptFile	<- list(
			YES=sprintf("%s/SURVIVEYES/tmp/GENES.txt",dataDir),
			NO=sprintf("%s/SURVIVENO/tmp/GENES.txt",dataDir)
		)
		# net ID-to-name mappings
		netInfo	<- list(
			YES=sprintf("%s/SURVIVEYES/tmp/NETWORKS.txt",dataDir),
			NO=sprintf("%s/SURVIVENO/tmp/NETWORKS.txt",dataDir)
			)
		# interaction nets
		netDir		<- list(
			YES=sprintf("%s/SURVIVEYES/tmp/INTERACTIONS",dataDir),
			NO=sprintf("%s/SURVIVENO/tmp/INTERACTIONS",dataDir)
		)
		# we are going to take union of FS pathways for each class so we need
		# the pathway scores for each class
		pTallyFile <- list(
			YES=sprintf("%s/SURVIVEYES/GM_results/SURVIVEYES_pathway_CV_score.txt",
						dataDir),
			NO=sprintf("%s/SURVIVENO/GM_results/SURVIVENO_pathway_CV_score.txt",
						dataDir)
		)
		
		# -----------------------------------------------------
		#if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
		#dir.create(outDir)

		# need pheno table for class assignment.
		pheno <- read.delim(clinical_file,sep="\t",h=T,as.is=T)
		colnames(pheno)[1] <- "ID"
		surv <- read.delim(survival_file,sep="\t",h=T,as.is=T)
		colnames(surv)[1:2] <- c("ID","STATUS_INT")
		survStr <- rep(NA,nrow(surv))
		survStr[surv$STATUS_INT<1] <- "SURVIVENO"
		survStr[surv$STATUS_INT>0] <- "SURVIVEYES"
		surv$STATUS <- survStr
		pheno <- merge(x=pheno,y=surv,by="ID")
		pheno$X <- NULL
		colnames(pheno)[which(colnames(pheno)=="STATUS")] <- "GROUP"
		
		# gene IDs should be the same for both classes
		tmp <- system(sprintf("diff %s %s", ptFile$YES,ptFile$NO),intern=TRUE)
		if (!is.null(attr(tmp,"status"))){
			cat("Gene IDs for all classes not identical! Assumption violated.\n")
			browser()
		}
		
		# pool best nets for both
		poolDir <- "./pool"
		if (file.exists(poolDir)) unlink(poolDir,recursive=TRUE)
		dir.create(poolDir)
		
		newNetIDs <- list()

		# pool feature selected nets from both groups
		outmat <- list() # num rows = num patients. # num cols = num unit
								     # features (e.g. one per gene). note this concatenates
										 # for all groups
		featTally <- list()
		for (gps in names(pTallyFile)) {
			cat(sprintf("Group %s\n", gps))
			pTally <- read.delim(pTallyFile[[gps]],sep="\t",h=T,as.is=T)
			pTally <- subset(pTally, pTally[,2]>=cutoff)[,1]
			#pTally <- sub("_cont|\\.profile","",pTally)
			if ("BRCA2" %in% pTally) pTally <- pTally[-which(pTally=="BRCA2")]
	
			featTally[[gps]] <- length(pTally)

			netInfo_cur <- read.delim(netInfo[[gps]],sep="\t",h=F,as.is=T)
			netInfo_cur[,2] <- sub("_cont|\\.profile","",netInfo_cur[,2])
			curNetIds <- matrix(NA,nrow=length(pTally),ncol=2)
			ctr <- 1

			# copy feature selected nets for this group
			for (cur in pTally) {
				if (verbose) cat(sprintf("%s\n",cur))

				if (any(grep(".profile", cur))) { # is an xpr matrix, store as-is
					if (verbose) cat("\t.profile net -> storing matrix\n")
					val <- read.delim(sprintf("%s/%s",dirname(netInfo[[gps]]),
						cur),sep="\t",h=F,as.is=T)
					rownames(val) <- val[,1]
					val <- val[,-1]

					if (length(outmat)>=1) {# is not first
						if (verbose) cat("not first, match ids\n")
						midx <- match(rownames(outmat[[1]]),rownames(val))
						val <- val[midx,,drop=FALSE]
						if (all.equal(rownames(val),rownames(outmat[[1]]))!=TRUE) {
							cat("indices don't match\n")
							browser()
					}}
					
					outmat[[cur]] <- as.matrix(val)
	
				} else if (any(grep("^MUT_",cur))) { # mutation-based net
					if (verbose) cat("\tMutation net -> storing as binary\n")
					cur <- sprintf("%s.txt",cur)
					tmp <- read.delim(sprintf("%s/%s",dirname(netInfo[[gps]]),
						cur),sep="\t",h=F,as.is=T)
				
					vec <- matrix(0,nrow=nrow(pheno),ncol=1)
					if (length(outmat)>=1) rownames(vec) <- rownames(outmat[[1]])
					else rownames(vec) <- pheno$ID
	
					hasmut <- union(tmp[,1],tmp[,2])
					vec[which(rownames(vec)%in% hasmut)] <- 1
					
					outmat[[cur]] <- vec
				} else {														 # clinical net
					if (verbose) cat("\tmust be clinical\n")
					if (cur=="Karnofsky_cont") cur2 <- "performance_score"
					else cur2 <- sub("_cont","",cur)
					val <- pheno[,cur2]
					if (class(val)!="numeric") {
						val <- as.integer(as.factor(val))
					}
					val <- as.matrix(val);
					rownames(val) <- pheno$ID
					
					if (length(outmat)>=1) {
						if (verbose) cat("\t\tnot first, match indices\n")
						midx <- match(rownames(outmat[[1]]),rownames(val))
						val <- val[midx,,drop=FALSE]
						if (all.equal(rownames(val),rownames(outmat[[1]]))!=TRUE) {
							cat("indices don't match\n")
							browser()
						}
					} 
					outmat[[cur]] <- val
					cur <- sprintf("%s.txt",cur)
				}
				idx <- which(netInfo_cur[,4] == cur)
				if (length(idx)<1) {
					cat(sprintf("%s: index not found!\n",cur))
					browser()
				}
		}
		}
		cat(sprintf("Compiled %i features\n", length(outmat)))
		print(featTally)

		allDat <- do.call("cbind",outmat)
		allDat <- t(allDat) # patients must becolumns for cosine
	
		# make net of cosine similarities.		
		cosSim <- lsa::cosine(x=na.omit(allDat))
		if (any(rowSums(is.na(cosSim))>0)) {
			cat("-> got NA\n")
			browser()
		}
		require(reshape2)
		cosSim <- melt(cosSim)
		cosSim[,1] <- as.character(cosSim[,1])
		cosSim[,2] <- as.character(cosSim[,2])
		idx <- which(cosSim[,1]==cosSim[,2]) # remove self edges
		cosSim <- cosSim[-idx,]

		netDx::sparsifyNet(cosSim, numPatients=ncol(allDat),outFile="tmp.txt",
				verbose=FALSE)
		tmp <- read.delim("tmp.txt",sep="\t",header=F,as.is=T)
		colnames(tmp) <- c("AliasA","AliasB","weight")
		write.table(tmp[,c(1,3,2)],file="tmp.txt",sep="\t",col=T,row=F,quote=F)
		
		aggNet <- tmp
		aggNet[,3] <- 1-aggNet[,3]  # convert similarity to dissimilarity

		cat("\n\n---------------------\n")
		cat(sprintf("Dijkstra distances: %s: %s :%s\n",
			curSet,aggFun,simMode))
		cat("---------------------\n")
		x <- compareShortestPath(aggNet, pheno,verbose=F)
		###print(lapply(x$avg, function(k) k[1]))
		
		pyes <- wilcox.test(x$all[["SURVIVEYES"]],
				x$all[["SURVIVEYES-SURVIVENO"]],
				alternative="less")$p.value
		pno <- wilcox.test(x$all[["SURVIVENO"]],
				x$all[["SURVIVEYES-SURVIVENO"]],
				 alternative="less")$p.value
		pyes2 <- wilcox.test(x$all[["SURVIVEYES"]],x$all[["overall"]],
				alternative="less")$p.value
		pno2 <- wilcox.test(x$all[["SURVIVENO"]],x$all[["overall"]],
				alternative="less")$p.value
		cat("WMW tests (one-sided, less)\n")
		cat(sprintf("YES: -> yes-no: p %1.2e ; -> ALL: p%1.2e\n", 
			pyes,pyes2))
		cat(sprintf("NO: -> yes-no: p %1.2e ; -> ALL: p%1.2e\n", 
			pno,pno2))

		curDijk[cur_i,] <- c(unlist(lapply(x$avg,function(k) k[1])),
			pyes,pyes2,pno,pno2)
		cur_i <- cur_i+1
		
		# plot density of shortest distances and compute i
		# statistics on significance
		if (!aggFun == "MAX") {
		pdf(sprintf("%s/%s_%s_%s_shortestDist.pdf",
						outDir, curSet,aggFun,simMode),width=8,height=6)
		tryCatch({
		par(bty='n',las=1,cex.axis=1.3)
		
		require(caroline)
		names(x$all) <- gsub("SURVIVE","",names(x$all))
		violins(x$all, deciles=FALSE,drawRect=TRUE, connect=c(), CImed=FALSE,
			ylab="Pairwise shortest path",ylim=c(0,0.8),
			col=c("orange","gray50","blanchedalmond","cyan"),las=1)
		#abline(h=0.5,lty=3,col='red')
		abline(h=median(x$all[[1]]),lty=1,lwd=2,col='black')
		# label violins with sample size
		for (k in 1:length(x$all)) {
				text(k,-0.01,
				sprintf("N=%s",prettyNum(length(x$all[[k]]),big.mark=',')),
				cex=1.3,font=3)
		}
		pvals <- matrix(nrow=4,ncol=4)
		my_y <-0.2
		for (i in 1:3) {
			for (j in (i+1):4) { 
					pvals[i,j] <- wilcox.test(x$all[[i]],x$all[[j]])$p.value
		
					# draw pvalue segments on violin plot
					segments(x0=i,x1=i,y0=my_y,y1=my_y+0.05,col='red')
					segments(x0=j,x1=j,y0=my_y,y1=my_y+0.05,col='red')
					segments(x0=i,x1=j,y0=my_y+0.05,y1=my_y+0.05,col='red')
		
					if (pvals[i,j] > 0.001) str <- sprintf("%1.3f",pvals[i,j])
					else str <- sprintf("p< %1.1e",pvals[i,j]) 
					text(x=(i+j)/2,y=my_y+0.03,str,labels=,cex=1.3,font=3)
					my_y <- my_y + 0.08
			}
		}
		},error=function(ex) {
			print(ex)
		},finally={
			dev.off()
		})
		}
		write.table(pheno,file=sprintf("%s/%s_%s_%s_PSN_pheno.txt",
			 outDir,curSet, simMode, aggFun),sep="\t",col=T,row=F,quote=F)

		# layout network in Cytoscape
		network.suid <- EasycyRest::createNetwork(
			nodes=pheno, nodeID_column="ID",edges=aggNet,
				netName=sprintf("%s_%s_%s",curSet,simMode,aggFun),
				collName=curSet
		)
		# spring-embedded layout on edge 'weight' column
		layout.url <- sprintf("%s/apply/layouts/kamada-kawai/%s?column=weight",
					base.url,network.suid, sep="/")
		response <- httr::GET(url=layout.url) 
		# apply style
		apply.style.url <- sprintf("%s/apply/styles/%s/%i",
					base.url,styleName,network.suid)
		response <- httr::GET(apply.style.url)
		
	}
	curDijk <- cbind(curDijk,method=sprintf("%s:%s", simMode,aggFun))
	dijk[[sprintf("%s_%s", aggFun,simMode)]] <- curDijk
}
}

dijk <- do.call("rbind",dijk)
write.table(dijk,file=sprintf("Dijkstra_PSN_%s.txt",dt),sep="\t",col=T,row=T,
	quote=F)
},error=function(ex){ 
	print(ex)
}, finally={
	sink(NULL)
})
		
