## ----eval=FALSE----------------------------------------------------------
## require(netDx)
## require(GenomicRanges)
## require(biomaRt) # for fetching gene coordinates
## 
## numCores	<- 8L
## set.seed(123)
## 
## outDir <- sprintf("%s/200129_threeWay",tempdir())
## if (file.exists(outDir)) unlink(outDir,recursive=TRUE);
## dir.create(outDir)
## 
## cat("* Setting up sample metadata\n")
## phenoFile <- sprintf("%s/extdata/AGP1_CNV.txt",
##                      path.package("netDx"))
## pheno   <- read.delim(phenoFile,sep="\t",header=T,as.is=T)
## colnames(pheno)[1] <- "ID"
## head(pheno)
## 
## cnv_GR    <- GRanges(pheno$seqnames,IRanges(pheno$start,pheno$end),
##                         ID=pheno$ID,LOCUS_NAMES=pheno$Gene_symbols)
## pheno <- pheno[!duplicated(pheno$ID),]
## 
## # create gene-range sets for pathways
## # first fetch pathway defintions and gene coordinates
## cat("* Create pathway range-sets\n")
## pathFile <- fetchPathwayDefinitions("February",2018,verbose=TRUE)
## pathwayList <- readPathways(pathFile)
## ensembl <- useMart("ENSEMBL_MART_ENSEMBL",
## 	dataset="hsapiens_gene_ensembl",
## 	host="may2009.archive.ensembl.org",
## 	path="/biomart/martservice",archive=FALSE)
## genes <- getBM(attributes=c("chromosome_name",
## 		"start_position",
## 		"end_position",
## 		"hgnc_symbol"),
## 	mart=ensembl)
## genes <- genes[which(genes[,4]!=""),]
## # create GRanges object and group by pathways
## gene_GR     <- GRanges(genes[,1],IRanges(genes[,2],genes[,3]),
##    name=genes[,4])
## path_GRList <- mapNamedRangesToSets(gene_GR,pathwayList)
## 
## # run the predictor
## predictClass	<- "case"
## nuildPredictor_sparseGenetic(pheno,predictClass,outDir,
## 			enrichLabels=TRUE,numPermsEnrich=20L,
## 			numCores=8L,featScoreMax=3L)
## 
## # evaluate performance
## plot(0,0,type="n",xlim=c(0,100),ylim=c(0,100),
##    las=1, xlab="FPR (%)", ylab="TPR (%)",bty='n',
## 	cex.axis=1.5)
## 	
## # pathway scores
## pathList <- list()
## 
## inFile <- sprintf("%s/RR_changeNetSum_stats_denEnrichedNets.txt",	
## 		outDir)
## dat	<- read.delim(inFile,sep="\t",header=TRUE,as.is=TRUE)
## points(dat$other_pct,dat$pred_pct,
## 	  col="red",type="o",pch=16,cex=0.5)
## 
## tmp <- data.frame(	
## 	score=dat$score,
## 	tp=dat$pred_ol,fp=dat$other_ol,
## 	# "-" that were correctly not called
## 	tn=dat$other_tot - dat$other_ol,
## 	# "+" that were not called
## 	fn=dat$pred_tot - dat$pred_ol)
## 
## stats <- netDx::perfCalc(tmp)
## tmp <- stats$stats
## cat(sprintf("PRAUC = %1.2f\n", stats$prauc))
## cat(sprintf("ROCAUC = %1.2f\n", stats$auc))
## 
## # now get pathway score
## ptmp <- read.delim(sprintf("%s/pathway_cumTally.txt",
##    outDir),sep="\t",h=T,as.is=T)


## ----eval=TRUE-----------------------------------------------------------
require(netDx)
require(GenomicRanges)
require(biomaRt) # for fetching gene coordinates

numCores	<- 8L 

outDir <- sprintf("%s/200129_threeWay",tempdir())
if (file.exists(outDir)) unlink(outDir,recursive=TRUE); dir.create(outDir)

cat("* Setting up sample metadata\n")
phenoFile <- sprintf("%s/extdata/AGP1_CNV.txt",path.package("netDx"))
pheno   <- read.delim(phenoFile,sep="\t",header=T,as.is=T)
colnames(pheno)[1] <- "ID"
head(pheno)

cnv_GR    <- GRanges(pheno$seqnames,IRanges(pheno$start,pheno$end),
                        ID=pheno$ID,LOCUS_NAMES=pheno$Gene_symbols)
pheno <- pheno[!duplicated(pheno$ID),]


## ----eval=TRUE-----------------------------------------------------------
cat("* Create pathway range-sets\n")
pathFile <- fetchPathwayDefinitions("February",2018,verbose=TRUE)
pathwayList <- readPathways(pathFile)

# get gene coordinates, use hg18
ensembl <- useMart("ENSEMBL_MART_ENSEMBL",
	dataset="hsapiens_gene_ensembl",
	host="may2009.archive.ensembl.org",
	path="/biomart/martservice",archive=FALSE)
genes <- getBM(attributes=c("chromosome_name",
		"start_position",
		"end_position",
		"hgnc_symbol"),
	mart=ensembl)
genes <- genes[which(genes[,4]!=""),]

gene_GR     <- GRanges(genes[,1],IRanges(genes[,2],genes[,3]),
   name=genes[,4])
path_GRList <- mapNamedRangesToSets(gene_GR,pathwayList)


## ----eval=TRUE-----------------------------------------------------------
predictClass	<- "case"
out <- buildPredictor_sparseGenetic(pheno, cnv_GR, predictClass,
                             path_GRList,outDir,
                             numSplits=3L, featScoreMax=3L,
                             enrichLabels=TRUE,numPermsEnrich=20L,
                             numCores=8L)


## ----eval=TRUE-----------------------------------------------------------
dat	<- out$performance_denEnrichedNets
plot(0,0,type="n",xlim=c(0,100),ylim=c(0,100),
	las=1, xlab="False Positive Rate (%)", 
	ylab="True Positive Rate (%)",
	bty='n',cex.axis=1.5,cex.lab=1.3,
	main="ROC curve - Patients in label-enriched pathways")
points(dat$other_pct,dat$pred_pct,
	  col="red",type="o",pch=16,cex=1.3,lwd=2)

tmp <- data.frame(	
	score=dat$score,
	tp=dat$pred_ol,fp=dat$other_ol,
	# tn: "-" that were correctly not called
	tn=dat$other_tot - dat$other_ol,
	# fn: "+" that were not called 
	fn=dat$pred_tot - dat$pred_ol) 

stats <- netDx::perfCalc(tmp)
tmp <- stats$stats
cat(sprintf("PRAUC = %1.2f\n", stats$prauc))
cat(sprintf("ROCAUC = %1.2f\n", stats$auc))

# examine pathway scores
print(out$cumulativeFeatscores)


