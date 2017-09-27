#' @param infoList (list) all information related to input data
#' to generate psn. includes:
#' 1. clinical_file: path to clinical_core.txt
#' 2. survival_file: path to binary_survival.txt
#' 3. outDir: where output files should be stored
#' 4. dataDir: root directory where input networks to be combined for psn
#' are stored
#' 5. ptFile: (list with YES and NO) path to GENES.txt file for each
#' 6. netInfo: (list with YES and NO) path to NETWORKS.txt file for each
#' category
#' 7. netDir: (list with YES and NO) path to INTERACTIONS/ dir
#' 8. netScoreFile (list with YES and NO) path to net score over 100 splits
#' 
getPSN <- function(infoList, consCutoff, consPctPass, topX,aggFun="MEAN",
	netMode="consensus",outDir=".") {

# --------------------------------------------------------------
# setup for network generation in Cytoscape
require(EasycyRest) # https://github.com/shraddhapai/EasycyRest
require(httr)
require(netDx)
require(ggplot2)

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

# --------------------------------------------
# work begins
for (nm in names(infoList)) {
	assign(nm, infoList[[nm]])
}

outPfx <- sprintf("%s_%s", setName, aggFun)
curDijk <- matrix(NA,nrow=1,ncol=8)
colnames(curDijk) <- c("SURVIVENO","SURVIVEYES","SURVIVEYES-SURVIVENO",
			"overall","pYES-YESNO",
			"pYES-overall","pNO-YESNO","pNO-overall")
# -----------------------------------------------------
#if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
#dir.create(outDir)
		
# gene IDs should be the same for both classes
if (length(ptFile)>=2) {
tmp <- system(sprintf("diff %s %s", ptFile$YES,ptFile$NO),intern=TRUE)
if (!is.null(attr(tmp,"status"))){
		cat("Gene IDs for all classes not identical! Assumption violated.\n")
		browser()
}
}

# pool best nets for both
poolDir <- sprintf("%s/pool",outDir)
if (file.exists(poolDir)) unlink(poolDir,recursive=TRUE)
dir.create(poolDir)
	newNetIDs <- list()

# pool feature selected nets from both groups
alreadyAdded <- c() 
for (gps in names(ptFile)) {
	cat(sprintf("Group %s\n", gps))
	
	netInfo_cur <- read.delim(netInfo[[gps]],sep="\t",h=F,as.is=T)
	netInfo_cur[,2] <- sub("_cont|\\.profile","",netInfo_cur[,2])
	
	if (setName %in% c("OV_oneClinNet","oneClinNet")) {
		pTally <- "clinical" # hard-coded for oneClinNet scenario
	} else if (setName %in% "oneRNANet") {
		pTally <- "rna"
	} else if (setName %in% "oneProtNet") {
		pTally <- "prot"
	} else {
		netScores	<- read.delim(netScoreFile[[gps]],sep="\t",h=T,as.is=T)
		netNames 	<- netScores[,1]
		netScores <- netScores[,-1]
		tmp <- rowSums(netScores >= consCutoff)
		idx <- which(tmp >= floor(consPctPass*ncol(netScores)))
		pTally <-  netNames[idx]
		pTally <- sub("_cont|\\.profile","",pTally)

		idx <- which(pTally %in% alreadyAdded)
		if (any(idx)) {
			cat(sprintf("Found a net already added before, removing: {%s}\n",
				paste(pTally[idx],collapse=",")))
			pTally <- pTally[-idx]
		}
		if (length(pTally)>=1) alreadyAdded <- c(alreadyAdded,pTally)
	}
	
	print(pTally)
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
		##cat(sprintf("\t%s -> %s\n", cur, tmp))
		ctr <- ctr+1
	}
	newNetIDs[[gps]] <- curNetIds
}

# write compiled net info file
newNetIDs <- do.call("rbind",newNetIDs)
netInfo_combinedF <- sprintf("%s/netInfo.txt", poolDir)
write.table(newNetIDs,file=netInfo_combinedF,sep="\t",
col=F,row=F,quote=F)

# now create integrated net using the best nets pooled 
# from all classes
cat(sprintf("%i networks with score >=%i\n",nrow(pTally),consCutoff))
cat("-----\n")

pdf(sprintf("%s/%s_edgedistr.pdf",outDir,outPfx))
aggNetFile <- netDx::writeWeightedNets(ptFile$YES,
	netInfo=netInfo_combinedF,
	poolDir,keepNets=newNetIDs[,2],outDir,
	filterEdgeWt=0,limitToTop=Inf,
	outFileName=sprintf("%s_PSN.txt",outPfx),
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


if (setName %in% "OV_oneClinNet") { ### OV pheno$GROUP results in 
			### creation of SURVIVENO-YES instead of YES-NO
			### which throws off the network layout method.
		  ### swap a couple rows to induce same ordering
			### as other tumour setes
			tmp <- pheno; tmp[1,] <- pheno[2,]; tmp[2,] <- pheno[1,];
			pheno <- tmp
		}	

write.table(pheno,file=sprintf("%s/%s_PSN_pheno.txt",
	outDir,outPfx),sep="\t",col=T,row=F,quote=F)

# get and plot Dijkstra distances for all group combinations
x <- compareShortestPath(aggNet, pheno,verbose=FALSE)
pdf(sprintf("%s/%s_DijkstraViolinplots.pdf", outDir,outPfx),width=10,height=4)
tryCatch({
	par(las=1,bty='n')
	dl <- data.frame(intType=rep(names(x$all),lapply(x$all,length)),
				dijk=unlist(x$all))
	plotList <- list()
	p <-ggplot(dl,aes(intType, dijk))
	p <- p + ylab("Pairwise Dijkstra distance\n(smaller is better)") 
	p <- p + xlab("Pair groups")
	p <- p + ggtitle(setName)
	p2 <- p+geom_violin(scale="width")+geom_boxplot(width=0.02) # + geom_jitter(width = 0.1,cex=0.3, alpha=0.5)
	print(p2)

	boxplot(x$all,pars=list(boxwex=0.4),cex.axis=0.8,
	xlab="Pair categories", ylab="Pairwise Dijkstra",
	title=sprintf("%s", outPfx))
},error=function(ex){
	print(ex)
},finally={
	dev.off()
})
for (k in 1:length(x$all)) {
	cat(sprintf("k=%i\n", k))

	idx <- which(colnames(curDijk) %in% names(x$all)[k])
	cat(sprintf("%s: idx=%i\n", names(x$all)[k], idx))
	curDijk[1,idx] <- median(x$all[[k]])
	cat(sprintf("%s: idx=%i, median = %1.2f\n", names(x$all)[k],idx,
		curDijk[1,idx]))
}

# compute Dijkstra p-values
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
cat(sprintf("YES: -> yes-no: p %1.2e ; -> ALL: p%1.2e\n", pyes,pyes2))
cat(sprintf("NO: -> yes-no: p %1.2e ; -> ALL: p%1.2e\n", pno,pno2))

curDijk[1,5:8] <-c(pyes,pyes2,pno,pno2)

# create a pruned network for visualization
aggNet_pruned <- netDx::pruneNetByStrongest(aggNet,pheno$ID, topX=topX)	
outFile <- sprintf("%s/%s_prunedNet_top%1.2f.txt",outDir,outPfx,topX)
write.table(aggNet_pruned,file=outFile,sep="\t",col=TRUE,row=FALSE,
	quote=FALSE)
# layout network in Cytoscape
network.suid <- EasycyRest::createNetwork(
	nodes=pheno, nodeID_column="ID",edges=aggNet_pruned,
	netName=sprintf("%s_%s_top%1.2f",setName,aggFun,topX),
	collName="KIRC"
)
# spring-embedded layout on edge 'weight' column
layout.url <- sprintf("%s/apply/layouts/kamada-kawai/%s?column=weight",
	base.url,network.suid, sep="/")
response <- httr::GET(url=layout.url) 
# apply style
apply.style.url <- sprintf("%s/apply/styles/%s/%i",
	base.url,styleName,network.suid)
response <- httr::GET(apply.style.url)

write.table(curDijk,
	file=sprintf("%s/%s_Dijkstra_PSN.txt",outDir,outPfx),
		sep="\t",col=T,row=T,quote=F)
}
