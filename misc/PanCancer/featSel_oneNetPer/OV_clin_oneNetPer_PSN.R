#' generate integrated PSN for OV - oneNetPer.
rm(list=ls())


# --------------------------------------------------------------
# Param for computing integrated PSN
consCutoff 		<-10  	# include nets with score >= this value
consPctPass		<- 1
dt <- format(Sys.Date(),"%y%m%d")

netMode <- "consensus" # consensus|bestAUC

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
	OV=sprintf("%s/2017_TCGA_OV/input/OV_clinical_core.txt",rootDir)
)
survList <- list(
	OV=sprintf("%s/2017_TCGA_OV/input/OV_binary_survival.txt",rootDir)
)

outDir <- sprintf("%s/2017_PanCancer_Survival/oneNetPer_FeatSel",rootDir)

cat("***** Consensus mode *****\n")
selIter <- list(
		OV=sprintf("%s/2017_TCGA_OV/output/OV_oneNetPer_170425", rootDir)
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

datSets <- "OV" # c("OV","OV","LUSC","GBM")
dijk <- list()
for (aggFun in c("MAX","MEAN")) {
	dijk[[aggFun]] <- list()
	for (simMode in "BinProp") { 		
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
			YES=sprintf("%s/tmp/GENES.txt",dataDir)
		)
		# net ID-to-name mappings
		netInfo	<- list(
			YES=sprintf("%s/tmp/NETWORKS.txt",dataDir)
			)
		# interaction nets
		netDir		<- list(
			YES=sprintf("%s/tmp/INTERACTIONS",dataDir)
		)
		
		# -----------------------------------------------------
		#if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
		#dir.create(outDir)
				
		# pool best nets for both
		poolDir <- "./pool"
		if (file.exists(poolDir)) unlink(poolDir,recursive=TRUE)
		dir.create(poolDir)
		
		newNetIDs <- list()

		# pool feature selected nets from both groups
		# running only for first group because both not needed
		for (gps in "YES") {
			cat(sprintf("Group %s\n", gps))
			
			netInfo_cur <- read.delim(netInfo[[gps]],sep="\t",h=F,as.is=T)
			netInfo_cur[,2] <- sub("_cont|\\.profile","",netInfo_cur[,2])
	
		if (netMode=="consensus") {
			cat("*** pTally hard-coded to clinical ***\n")
			pTally <- "clinical"
			print(pTally)
		}
			curNetIds <- matrix(NA,nrow=length(pTally),ncol=2)
			ctr <- 1

			# copy feature selected nets for this group
			for (cur in pTally) {
				idx <- which(netInfo_cur[,2] == cur)
				if (length(idx)<1) {
					cat(sprintf("%s: index not found!\n",cur))
					browser()
				}
				netID <- sprintf("1.%s.txt",netInfo_cur[idx,1])
				tmp <- sprintf("%s.%s", gps,netID)
				file.copy(from=sprintf("%s/%s", netDir[[gps]], netID), 
							to=sprintf("%s/1.%s",poolDir,tmp))
					
				curNetIds[ctr,] <- c(sub(".txt","",tmp), 
									 paste(gps,netInfo_cur[idx,2],sep=".")) 
				ctr <- ctr+1
			}
			newNetIDs[[gps]] <- curNetIds
		}

		# write compiled net info file
		newNetIDs <- do.call("rbind",newNetIDs)
		if (simMode == "BinProp") {
			isBinary <- rep(0,nrow(newNetIDs))
			newNetIDs <- cbind(newNetIDs, isBinary=isBinary)
		}
		netInfo_combinedF <- sprintf("%s/netInfo.txt", poolDir)
		write.table(newNetIDs,file=netInfo_combinedF,sep="\t",
			col=F,row=F,quote=F)
		
		# now create integrated net using the best nets pooled 
		# from all classes
		cat(sprintf("%i networks with score >=%i\n",nrow(pTally),consCutoff))
		cat("-----\n")
pdf("test.pdf")
		aggNetFile <- netDx::writeWeightedNets(ptFile$YES,
						netInfo=netInfo_combinedF,
						poolDir,keepNets=newNetIDs[,2],outDir,
						filterEdgeWt=0,limitToTop=Inf,
						outFileName=sprintf("%s_%s_%s_PSN.txt", 
							curSet,simMode,aggFun),
						plotEdgeDensity=TRUE,
						writeAggNet=aggFun,verbose=FALSE)
dev.off()
		aggNet<- read.delim(aggNetFile,sep="\t",h=T,as.is=T)[,1:3]
		colnames(aggNet) <- c("AliasA","AliasB","weight")
		aggNet[,3] <- 1-aggNet[,3]  # convert similarity to dissimilarity
		
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
		cat("\n\n---------------------\n")
		cat(sprintf("Dijkstra distances: %s: %s :%s\n",
			curSet,aggFun,simMode))
		cat("---------------------\n")
		x <- compareShortestPath(aggNet, pheno,verbose=FALSE)
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
		write.table(pheno,file=sprintf("%s/%s_%s_%s_%s_PSN_pheno.txt",
			 outDir,curSet, simMode, aggFun,netMode),sep="\t",col=T,row=F,quote=F)

		# create a pruned network for visualization
		aggNetFile <- netDx::writeWeightedNets(ptFile$YES,
						netInfo=netInfo_combinedF,
						poolDir,keepNets=newNetIDs[,2],outDir,
						filterEdgeWt=0,limitToTop=25,
						outFileName=sprintf("%s_%s_%s_PSNpruned.txt", 
							curSet,simMode,aggFun),
						writeAggNet=aggFun,verbose=FALSE)
		aggNet_pruned<- read.delim(aggNetFile,sep="\t",h=T,as.is=T)[,1:3]
		colnames(aggNet_pruned) <- c("AliasA","AliasB","weight")
		# convert similarity to dissimilarity
		aggNet_pruned[,3] <- 1-aggNet_pruned[,3]  

		# layout network in Cytoscape
		network.suid <- EasycyRest::createNetwork(
			nodes=pheno, nodeID_column="ID",edges=aggNet_pruned,
				netName=sprintf("%s_clinOneNet_%s",curSet,aggFun),
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
write.table(dijk,
	file=sprintf("%s/%s_Dijkstra_PSN_%s_%s.txt",
	outDir,curSet,aggFun,dt),
	sep="\t",col=T,row=T,quote=F)

},error=function(ex){ 
	print(ex)
}, finally={
	sink(NULL)
})
		
