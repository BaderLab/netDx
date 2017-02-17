#' generate integrated PSN for TCGA breast cancer data
rm(list=ls())

# to create pheno table
clinical_file <- "/mnt/data2/BaderLab/PanCancer_GBM/input/GBM_clinical_core.txt"            
survival_file <- "/mnt/data2/BaderLab/PanCancer_GBM/input/GBM_binary_survival.txt"            

# directory with output to process
dataDir <- "/mnt/data2/BaderLab/PanCancer_GBM/output/featSel_170209/rng4"
#dataDir <- "/mnt/data2/BaderLab/PanCancer_GBM/output/featSel_noMut_170213/rng7"
# where the integrated net will be written
outDir		<- "."

cutoff 		<-10  # include nets with score >= this value
corrThresh 	<-0.7 # exclude edges with similarity lower than this threshold
aggFun 		<-"MAX" # aggregate edges between any given pair of patients

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
require(netDx)

# gene IDs should be the same for both classes
tmp <- system(sprintf("diff %s %s", ptFile$YES,ptFile$NO),intern=TRUE)
if (length(tmp)>0) {
	cat("Gene IDs for all classes not identical! Assumption violated.\n")
	browser()
}

# pool best nets for both
poolDir <- "./pool"
if (file.exists(poolDir)) unlink(poolDir,recursive=TRUE)
dir.create(poolDir)

newNetIDs <- list()
for (gps in names(pTallyFile)[2]) {
	cat(sprintf("Group %s\n", gps))
	pTally <- read.delim(pTallyFile[[gps]],sep="\t",h=T,as.is=T)
	pTally <- subset(pTally, pTally[,2]>=cutoff)[,1]
	pTally <- sub("_cont|\\.profile","",pTally)
	netInfo_cur <- read.delim(netInfo[[gps]],sep="\t",h=F,as.is=T)
	netInfo_cur[,2] <- sub("_cont","",netInfo_cur[,2])
	# copy feature selected nets for this group
	curNetIds <- matrix(NA,nrow=length(pTally),ncol=2)
	ctr <- 1
	for (cur in pTally) {
		idx <- which(netInfo_cur[,2] == cur)
		if (length(idx)<1) cat(sprintf("%s: index not found!\n",cur))
		netID <- sprintf("1.%s.txt",netInfo_cur[idx,1])
		tmp <- sprintf("%s.%s", gps,netID)
		file.copy(from=sprintf("%s/%s", netDir[[gps]], netID), 
					to=sprintf("%s/1.%s",poolDir,tmp))
			
		curNetIds[ctr,] <- c(sub(".txt","",tmp), 
							 paste(gps,netInfo_cur[idx,2],sep=".")) 
		cat(sprintf("\t%s -> %s\n", cur, tmp))
		ctr <- ctr+1
	}
	newNetIDs[[gps]] <- curNetIds
}
# write compiled net info file
newNetIDs <- do.call("rbind",newNetIDs)
netInfo_combinedF <- sprintf("%s/netInfo.txt", poolDir)
write.table(newNetIDs,file=netInfo_combinedF,sep="\t",col=F,row=F,quote=F)

# now create integrated net using the best nets pooled from all classes
cat(sprintf("%i networks with score >=%i\n",nrow(pTally),cutoff))
cat("-----\n")
aggNetFile <- netDx::writeWeightedNets(ptFile$YES,
				netInfo=netInfo_combinedF,
				poolDir,keepNets=newNetIDs[,2],outDir,
				filterEdgeWt=corrThresh,limitToTop=25,
				writeAggNet=aggFun,verbose=TRUE)
aggNet<- read.delim(aggNetFile,sep="\t",h=T,as.is=T)[,1:3]
colnames(aggNet)[3] <- "weight"

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
pheno$gender <- ifelse(pheno$gender=="FEMALE",1, 0)

colnames(pheno)[which(colnames(pheno)=="STATUS")] <- "GROUP"
aggNet[,3] <- 1-aggNet[,3]  # convert similarity to dissimilarity
x <- compareShortestPath(aggNet, pheno)

# plot density of shortest distances and compute statistics on significance
pdf(sprintf("%s/GBM_shortestDist.pdf",outDir),width=8,height=6)
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

write.table(pheno,file=sprintf("%s/PSN_pheno.txt",outDir),sep="\t",
			col=T,row=F,quote=F)


