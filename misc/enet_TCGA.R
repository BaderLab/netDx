# run TCGA example by prediction using elastic net
# Note that this example runs out of stack memory when run in RStudio.
# On her purple Mac, SP finds that it runs when you start R from command
# line, using the --max-ppsize=500000 option
# $ R --max-ppsize=500000

# caret is the wrapper package with functions to streamline
# model building
require(caret)

# load the data
require(netDx.examples)
data(TCGA_BRCA)

# convert into binary problem
pheno$STATUS[which(!pheno$STATUS %in% "LumA")] <- "other"
table(pheno$STATUS)
# train() expects patients in rows,genes in columns
xpr <- as.data.frame(t(xpr)) 
xpr$STATUS <- as.factor(pheno$STATUS)
                
print(all.equal(rownames(xpr),pheno$ID))

# partition into training and test
cat("* Creating training/test set\n")
inTrain <- createDataPartition(
    y=xpr$STATUS,p=0.75, list=FALSE)
training <- xpr[inTrain,]
testing <- xpr[-inTrain,]

### train the elastic net
# first set up the resampling
# this trainControl is taken from the kaggle example here:
# http://www.r-bloggers.com/kaggle-competition-walkthrough-fitting-a-model/
trControl <- trainControl(method="repeatedCV",
    number=10,repeats=3,classProbs=TRUE,
    summaryFunction=twoClassSummary)
# train elastic net and include parameter tuning
cat("* Running elastic net model builder ...\n")
enet <- train(STATUS~., data=training,
    method="glmnet",metric="ROC",trControl=trControl)
cat("done\n")

# assess model performance
## we have to exclude 43 samples because of missing data
numsamp <- nrow(testing)
testing <- na.omit(testing)
cat(sprintf("Removing missing data cuts N down from %i to %i\n",
	numsamp,nrow(testing)))
predClass <- predict(enet,newdata=testing)
confusionMatrix(predClass, testing$STATUS)

