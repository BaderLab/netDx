#' generate integrated PSN for PanCancer dataset
#' 

rm(list=ls())
require(EasycyRest) # https://github.com/shraddhapai/EasycyRest
require(httr)
require(netDx)

# --------------------------------------------------------------
# Param for computing integrated PSN
corrThresh 	<-0.7 	# exclude edges with similarity < this threshold
dt 					<- format(Sys.Date(),"%y%m%d")
aggFun		<- "MAX"

# --------------------------------------------------------------
# setup for network generation in Cytoscape

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
	LUSC=sprintf("%s/2017_TCGA_LUSC/input/LUSC_clinical_core.txt",rootDir)
)
survList <- list(
	LUSC=sprintf("%s/2017_TCGA_LUSC/input/LUSC_binary_survival.txt",rootDir)
)

outDir <- sprintf("%s/2017_PanCancer_Survival/integratedPSN",rootDir)

selIter <- list(
		LUSC=sprintf("%s/2017_TCGA_LUSC/output/ownTrain_170205",rootDir)
)

# which nets to combine for a given cancer type. selected based on 
# performance evaluate independently.
pTallySet <- list(
	LUSC=c("rppa")
)


# --------------------------------------------------------------
# Work begins

logFile <- sprintf("%s/getPSN_%s.log",outDir,dt)
sink(logFile,split=TRUE)

tryCatch({

datSets <- c("OV","KIRC","LUSC","GBM")
dijk <- list()
		curDijk <- matrix(NA,nrow=length(datSets),ncol=8)
		rownames(curDijk) <- datSets
		colnames(curDijk) <- c("NO","YES","YES-NO","overall","pYES-YESNO",
			"pYES-overall","pNO-YESNO","pNO-overall")

		cur_i <- 1

		for (curSet in "LUSC") {		

		cat(sprintf("%s\n", curSet))
		
		clinical_file <- clinList[[curSet]]
		survival_file <- survList[[curSet]]
		dataDir <- selIter[[curSet]]
		ptFile <- sprintf("%s/tmp/GENES.txt",dataDir)
		netInfo <- sprintf("%s/tmp/NETWORKS.txt",dataDir)
		netDir	<- sprintf("%s/tmp/INTERACTIONS",dataDir)

		# pool best nets for both
		poolDir <- sprintf("%s/pool",outDir)
		if (file.exists(poolDir)) unlink(poolDir,recursive=TRUE)
		dir.create(poolDir)
		
		newNetIDs <- list()

		# pool feature selected nets from both groups
		netInfo_cur <- read.delim(netInfo,sep="\t",h=F,as.is=T)
		netInfo_cur[,2] <- sub("_cont|\\.profile","",netInfo_cur[,2])
		pTally <- pTallySet[[curSet]]
	
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
				file.copy(from=sprintf("%s/%s", netDir, netID), 
							to=sprintf("%s/1.%s",poolDir,netID))
				curNetIds[ctr,] <- sub(".txt","",netID)
				ctr <- ctr+1
			}
			newNetIDs <- curNetIds

		# write compiled net info file
		netInfo_combinedF <- sprintf("%s/netInfo.txt", poolDir)
		write.table(newNetIDs,file=netInfo_combinedF,sep="\t",
			col=F,row=F,quote=F)
		
		# now create integrated net using the best nets pooled 
		# from all classes
		aggNetFile <- netDx::writeWeightedNets(ptFile,
						netInfo=netInfo_combinedF,
						poolDir,keepNets=newNetIDs[,2],outDir,
						filterEdgeWt=corrThresh,limitToTop=25,
						outFileName=sprintf("%s_noFeatSel_PSN.txt",curSet),
						writeAggNet=aggFun,verbose=FALSE)
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
		cat(sprintf("Dijkstra distances: %s: %s\n",
			curSet,aggFun))
		cat("---------------------\n")
		x <- compareShortestPath(aggNet, pheno,verbose=F)
		
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
		
			write.table(pheno,file=sprintf("%s/%s_%s_PSN_pheno.txt",
			 outDir,curSet, aggFun),sep="\t",col=T,row=F,quote=F)

		# layout network in Cytoscape
		network.suid <- EasycyRest::createNetwork(
			nodes=pheno, nodeID_column="ID",edges=aggNet,
				netName=sprintf("%s_%s",curSet,aggFun),
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

#write.table(dijk,
#	file=sprintf("%s/Dijkstra_PSN_%s_%s_%s_%s.txt",
#	outDir,simMode,aggFun,netMode,dt),
#	sep="\t",col=T,row=T,quote=F)

},error=function(ex){ 
	print(ex)
}, finally={
	sink(NULL)
})
		
