# run TCGA example by prediction using elastic net
# Note that this example runs out of stack memory when run in RStudio.
# On her purple Mac, SP finds that it runs when you start R from command
# line, using the --max-ppsize option
# $ R --max-ppsize=5000000 (six zeroes)

# Instructions for use of methods were obtained from the following pages:
# Elastic net: http://www.r-bloggers.com/kaggle-competition-walkthrough-fitting-a-model/
# Random forest: http://bigcomputing.blogspot.ca/2014/10/an-example-of-using-random-forest-in.html
rm(list=ls())


outDir <- "~/tmp/netDX/results"
method2use <- "RandomForest" ###ElasticNet | RandomForest
rmMissing  <- FALSE
rngSeed <- 102
set.seed(rngSeed) # make reproducible

if (!file.exists(outDir)) dir.create(outDir,recursive=TRUE)

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
require(netDx.examples)
data(TCGA_BRCA)

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
mega_TT <- rep("TEST",nrow(pheno))
mega_TT[inTrain] <- "TRAIN" ## pre-assigned
pheno$TT_STATUS <- mega_TT
print(table(pheno[,c("STATUS","TT_STATUS")]))

if (rmMissing) {
	nr <- nrow(xpr)
	cat("*** Excluding missing data\n")
	xpr <- na.omit(xpr)
	cat(sprintf("%i of %i genes excluded\n",nr-nrow(xpr), nr))
}
# train() expects patients in rows,genes in columns
xpr <- as.data.frame(t(xpr)) 
xpr$STATUS <- as.factor(pheno$STATUS)
                
print(all.equal(rownames(xpr),pheno$ID))

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
	mod <- train(STATUS~., data=training,
	    method="glmnet",metric="ROC",trControl=trControl)
	cat("done\n")
	
## random forest
} else if (method2use=="RandomForest") {

	cat("* Training random forest\n")
	mod <- train(STATUS~., data=training, 
		method="rf",
		trControl=trainControl(method="repeatedCV",number=10,repeats=3),
		prox=TRUE,allowParallel=TRUE);
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
	stopCluster(cl)
})


