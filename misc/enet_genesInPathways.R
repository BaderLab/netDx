# run TCGA example by prediction using elastic net
# Note that this example runs out of stack memory when run in RStudio.
# On her purple Mac, SP finds that it runs when you start R from command
# line, using the --max-ppsize option
# $ R --max-ppsize=5000000 (six zeroes)
#
# VARIATION 1 -- limit to genes in pathways

# Instructions for use of methods were obtained from the following pages:
# Elastic net: http://www.r-bloggers.com/kaggle-competition-walkthrough-fitting-a-model/
# Random forest: http://bigcomputing.blogspot.ca/2014/10/an-example-of-using-random-forest-in.html
args <- commandArgs(TRUE)
method2use <- args[1] ###ElasticNet | RandomForest

#outDir <- "/Users/shraddhapai/Google Drive/PatientNetworks/Papers/netDX/results/genes_in_pathways"
outDir <- "~/tmp/netDX/results/genes_in_pathways"
if (!file.exists(outDir)) dir.create(outDir)

rmMissing  <- TRUE
rngSeed <- 102
set.seed(rngSeed) # make reproducible

require(foreach)
require(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl)

# ---------------------------------------------------
# setup
dt <- format(Sys.Date(),"%y%m%d")
fPrefix <- sprintf("%s/%s_rmMiss%s_%s", outDir, method2use, rmMissing,dt)
logFile <- sprintf("%s.log",fPrefix)

# caret is the wrapper package with functions to streamline
# model building
require(caret)

# load the data
require(netDx)
require(netDx.examples)
data(TCGA_BRCA)

# get pathways 
pathFile 	<- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
    path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)
pathGenes	<- unlist(pathwayList)

#### create partition once.
###inTrain <- createDataPartition(
###    y=pheno$STATUS,p=0.67, list=FALSE)
###write.table(inTrain,file=sprintf("%s/trainID.txt",outDir),sep="\t",
###			col=F,row=F,quote=F)
###return(T)

sink(logFile,split=TRUE)
tryCatch({
print(Sys.time())

# read in preassigned train/test split
inTrain <- read.delim(sprintf("%s/trainID.txt",outDir),sep="\t",h=F,as.is=T)
inTrain <- inTrain[,1]

# convert into binary problem
pheno$STATUS[which(!pheno$STATUS %in% "LumA")] <- "other"
mega_TT 			<- rep("TEST",nrow(pheno))
mega_TT[inTrain] 	<- "TRAIN" ## pre-assigned
pheno$TT_STATUS 	<- mega_TT
print(table(pheno[,c("STATUS","TT_STATUS")]))

if (rmMissing) {
	nr <- nrow(xpr)
	cat("*** Excluding missing data\n")
	xpr <- na.omit(xpr)
	cat(sprintf("%i of %i genes excluded\n",nr-nrow(xpr), nr))
}

cat("**** Limiting to genes in pathways: ")
tmp	<- nrow(xpr)
xpr <- xpr[which(rownames(xpr) %in% pathGenes),]
cat(sprintf("gene xpr: genes in pathways %i --> %i after filter\n", 
			tmp,nrow(xpr)))

# add CNV data 
names(pathGenes) <- NULL
pathGenes <- unique(pathGenes)
data(genes)
genes <- subset(genes, name2 %in% pathGenes)
# combine coordinates for multiple gene isoforms on the same chromosome
# so that the extent is as wide as possible
x <- aggregate(genes$txStart,
			   by=list(gene=genes$name2,chrom=genes$chrom),
			   FUN=min)
y <- aggregate(genes$txEnd, 
			   by=list(gene=genes$name2,chrom=genes$chrom), 
			   FUN=max)
z <- merge(x,y,by=c("gene","chrom"))
cat(sprintf("After merging, %i genes -> %i merged\n", nrow(genes),nrow(z)))
genes <- z; rm(x,y,z)

gene_GR     <- GRanges(genes$chrom,
   IRanges(genes$x.x,genes$x.y),name=genes$gene)
# limit genes to those on chroms for which cnvs exist
seqlevels(gene_GR,force=TRUE) <- seqlevels(cnv_GR)
cat(sprintf("after dropping extra levels, have %i genes\n",length(gene_GR)))

# add cnv data
# create matrix with binary variables, one per pathway
# xpr3[i,j]=1 if patient j has cnv in pathway i; else 0
require(foreach)
require(parallel)
require(bigmemory)
cl <- makeCluster(6) # add outfile="" to print all worker output to stdout
registerDoParallel(cl)
bkFile <- sprintf("%s/tmp.bk",outDir)
if (file.exists(bkFile)) unlink(bkFile)
xpr3 <- big.matrix(NA,nrow=length(gene_GR),ncol=ncol(xpr),
				   type="integer",backingfile="tmp.bk",
				   backingpath=outDir,descriptorfile="tmp.desc",
				   dimnames=list(NULL,colnames(xpr)))
t0 <- Sys.time()
x <- foreach (k=1:length(gene_GR), .packages=c("GenomicRanges")) %dopar% {
	seqlevels(gene_GR[k], force=TRUE) <- seqlevels(cnv_GR)
	ol <- findOverlaps(cnv_GR, gene_GR[k])
	ol <- ol@queryHits
	if (length(ol)>0) {
		samps <- which(colnames(xpr) %in% unique(cnv_GR$ID[ol]))
		m <- bigmemory::attach.big.matrix(sprintf("%s/tmp.desc",outDir))
		m[k,] <- 0L
		m[k,samps] <- 1L
	}
}
print(Sys.time()-t0)
stopCluster(cl)
registerDoSEQ()
xpr3 <- as.matrix(xpr3)
xpr3 <- na.omit(xpr3)
cat(sprintf("gene-level CNV: Added %i terms\n", nrow(xpr3)))

# train() expects patients in rows,genes in columns
xpr <- t(xpr)
xpr3 <- t(xpr3)
if (all.equal(rownames(xpr),rownames(xpr3))!=TRUE) {
	cat("rownames don't match\n")
	browser()
}
xpr <- cbind(xpr,xpr3)
xpr <- as.data.frame(xpr)

if (all.equal(rownames(xpr),pheno$ID)!=TRUE) {
	cat("pheno, xpr rownames don't match");
	browser()
}
xpr$STATUS <- as.factor(pheno$STATUS)
                
cat("Final matrix submitted for ML:\n")
print(dim(xpr))

# ---------------------------------------------------
# model fitting

# partition into training and test
cat("* Creating training/test set\n")
training <- xpr[inTrain,]
testing <- xpr[-inTrain,]

## elastic net
if (method2use == "ElasticNet") {
	# first set up the resampling
	trControl <- trainControl(method="repeatedCV",
	    number=10,repeats=3,classProbs=TRUE,
	    summaryFunction=twoClassSummary)
	# train elastic net and include parameter tuning
	cat("* Running elastic net model builder ...\n")
	t0 <- Sys.time()
	mod <- train(STATUS~., data=training,
	    method="glmnet",metric="ROC",trControl=trControl)
	cat("done\n")
	print(Sys.time()-t0)
	
## random forest
} else if (method2use=="RandomForest") {

	cat("* Training random forest\n")
	t0 <- Sys.time()
	mod <- train(STATUS~., data=training, 
		method="rf",
		trControl=trainControl(method="repeatedCV",number=10,repeats=3),
		prox=TRUE,allowParallel=TRUE);
	print(Sys.time()-t0)
}

# -----------------------------------------------------
# predictor evaluation

# remove missing samples, otherwise the confusionMatrix() call
# throws an error. predict() quietly removes missing values, 
# so that confusionMatrix() is now comparing a pared-down version
# of test samples with the full version
nr <- nrow(testing)
testing_nona <- na.omit(testing)
df <- nr-nrow(testing_nona)
cat(sprintf("%i of %i samples ignored because of NA (%i%%)\n", 
			df, nr, round((df/nr)*100)))
testing <- testing_nona

predClass 	<- predict(mod,newdata=testing)
confmat		<- confusionMatrix(predClass, testing$STATUS)
# get variable importance
x <- varImp(mod,scale=FALSE)
y <- as.numeric(x$importance[,1]); names(y) <- rownames(x$importance)
y <- y[order(y,decreasing=TRUE)]

print(confmat)

dt <- format(Sys.Date(),"%y%m%d")
out <- list(method=method2use, rmMissing=rmMissing,
			mod=mod, predClass=predClass,testing=testing,
			inTrain=inTrain,RNGseed=rngSeed,
			confmat=confmat)
save(out,file=sprintf("%s.Rdata",fPrefix))
#plot(y,xlim=c(0,50),type='n',xlab="Rank",ylab="Coeff",
#	 main="TCGA BRCA: Gene weights from elastic net",las=1,
#	 bty='n')
#text(x=1:30,y=y[1:30],names(y)[1:30],font=3)
},error=function(ex){
	print(ex)
},finally={
	print(Sys.time())
	sink(NULL)
})


