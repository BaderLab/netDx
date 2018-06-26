#' look at YES/NO separation when only using individual clin nets
rm(list=ls())
require(igraph)

rootDir <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/"
# directory with inidividual sim nets
inDir <- sprintf("%s/2017_TCGA_KIRC/output/KIRC_clinNets_170430/rng1/SURVIVEYES/tmp",
	rootDir)

outDir <- sprintf("%s/2017_PanCancer_Survival/clinNets_170430",rootDir)
phenoFile <- sprintf("%s/2017_TCGA_KIRC/input/KIRC_binary_survival.txt",rootDir)

# pruned net
#cytoNetFile <- sprintf("%s/KIRC_BinProp_MEAN_PSNpruned.txt",
#	outDir)
cytoNetFile <- sprintf("%s/PSN_current/KIRC_BinProp_MEAN_PSN.txt",outDir)

# -----------------------------------------------------------
# helper 
.getAvgD <- function(mat) {
		tmp <- mat[upper.tri(mat,diag=FALSE)]
		idx <- which(is.infinite(tmp))
		if (any(idx)) tmp <- tmp[-idx]
		##if (verbose) cat(sprintf("\tN=%i distances\n", length(tmp)))

		c(mean(tmp,na.rm=TRUE), sd(tmp,na.rm=TRUE),length(tmp))
	}


# -----------------------------------------------------------
nets <- list(age="age_cont.txt",
	stage="stage_cont.txt",grade="grade_cont.txt")

pheno <- read.delim(phenoFile,sep="\t",h=T,as.is=T)
surv <- rep(NA,nrow(pheno))
surv[which(pheno[,2]==0)] <- "SURVIVENO"
surv[which(pheno[,2]==1)] <- "SURVIVEYES"
pheno$GROUP <- surv
colnames(pheno) <- c("ID","tmp","GROUP")

# plot distribution 
yes <- pheno$ID[which(pheno$GROUP=="SURVIVEYES")]
no <- pheno$ID[which(pheno$GROUP=="SURVIVENO")]
cytoNet <- read.delim(cytoNetFile,sep="\t",h=T,as.is=T)[,1:3]
cytoNet$value <- 1-cytoNet$value # convert to dissimilarity

.showDijk <- function(net,ttl="dijk") {
# get Dijkstra
	x <- netDx::compareShortestPath(net[,1:3],pheno,verbose=F);
	blah <- x$all
	y <- wilcox.test(blah$SURVIVEYES, 
		blah[["SURVIVENO-SURVIVEYES"]],alternative="less")
	z <- wilcox.test(blah$SURVIVENO, 
		blah[["SURVIVENO-SURVIVEYES"]],alternative="less")
	cat(sprintf("YY<YN: p < %1.2e ; NN < YN: p < %1.2e \n", 
		y$p.value,z$p.value))
	
	dl <- data.frame(intType=rep(names(blah),sapply(blah,length)), 
		dijk=unlist(blah))
	
	require(ggplot2)
	source("../../multiplot.R") # plot multiple ggplot panels on one fig.
	plotList <- list()
	p <-ggplot(dl,aes(intType, dijk))
	p <- p + ylab("Pairwise Dijkstra distance\n(smaller is better)") 
	p <- p + xlab("Pair groups")
	p <- p + ggtitle(ttl)
	p2 <- p+geom_boxplot() + geom_jitter(width = 0.1,cex=0.3, alpha=0.5)
	print(p2)
}

# ----------------------
# prune for top X edges and plot in Cytoscape
colnames(cytoNet)[3] <- "weight"
source("pruneByStrongest.R")
cy2 <- pruneByStrongest(cytoNet,pheno$ID,topX=0.2)
.showDijk(cy2,ttl="pruned, greatest 1% distances")
colnames(cy2)[1:3] <- c("AliasA","AliasB","weight")
cy2[,1] <- as.character(cy2[,1])
cy2[,2] <- as.character(cy2[,2])

source("makePSNstyle.R") # set up cytoscape
	# layout network in Cytoscape
network.suid <- EasycyRest::createNetwork(
	nodes=pheno, nodeID_column="ID",edges=cy2,
	netName="pruned",collName="KIRC2"
)
# spring-embedded layout on edge 'weight' column
layout.url <- sprintf("%s/apply/layouts/kamada-kawai/%s?column=weight",
base.url,network.suid, sep="/")
response <- httr::GET(url=layout.url) 
# apply style
apply.style.url <- sprintf("%s/apply/styles/%s/%i",
base.url,styleName,network.suid)
response <- httr::GET(apply.style.url)

# ----------------------
browser()


# run dijkstra for each net in turn
outmat <- matrix(NA,nrow=length(nets),ncol=6)
colnames(outmat) <- c("no mean","yes mean", "yes-no mean", "overall mean",
		"no_lt_yesno","yes_lt_yesno")
rownames(outmat) <- names(nets)
ctr <- 1
for (nm in names(nets)) {
	print(nm)
	net <- read.delim(sprintf("%s/%s", inDir,nets[[nm]]),sep="\t",h=F,as.is=T)
	net[,3] <- 1-net[,3] # make dissimilarity
	colnames(net) <- c("SOURCE","TARGET","WEIGHTS")
	write.table(net,file=sprintf("%s/%s_dissimilarity.txt",outDir,nm),sep="\t",
		col=T,row=F,quote=F)

	out <- netDx::compareShortestPath(net,pheno)

	outmat[ctr,1:4] <- unlist(lapply(out$avg, function(x) x[1]))
	outmat[ctr,5] <- wilcox.test(out$all[["SURVIVENO"]],
				out$all[["SURVIVENO-SURVIVEYES"]],"less")$p.value
	outmat[ctr,6] <- wilcox.test(out$all[["SURVIVEYES"]],
				out$all[["SURVIVENO-SURVIVEYES"]],"less")$p.value
	ctr <- ctr+1
}

write.table(outmat,file=sprintf("%s/clinNets_indivDijk.txt",outDir),sep="\t",
	col=T,row=T,quote=F)
