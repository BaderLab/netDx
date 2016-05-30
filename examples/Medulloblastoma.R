#' 4-way classifier built with GeneMANIA using the medulloblastoma dataset
#' from Northcott et al. (2012) Acta Neuropathologica. 
#' The paper above identified 4 subgroups of patients,
#' each predictable by the expression of 5-6 signature genes. 
#' Here we predict the class of a new patient using GeneMANIA.
#' Conceptually, we build 4 classifiers, one per subtype. We get 
#' get patient ranking for each of the classifier and assign the class  
#' based on which predictor gives the patient the highest rank.
#'
rm(list=ls())

# change this to a directory to which you have write access
outDir <-  "~/tmp/MB"
dir.create(outDir)

numCores    <- 2L	# number of cores for parallel processing
pctTrain    <- 0.7	# fraction of samples to use for feature selection

require(netDx)
require(netDx.examples)

# Load the Medulloblastoma dataset
data(MBlastoma)

# subtypes and genes predictive of these. From Table 1 of PMC3306784
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3306784/table/Tab1/
groupSig <- list(
    WNT=c("WIF1","TNC","GAD","DKK2","EMX2"),
    SHH=c("PDLIM3","EYA1","HHIP","ATOH1","SFRP1"),
    Group3=c("IMPG2","GABRA5","EGFL11","NRL","MAB21L2","NPR3"),
    Group4=c("KCNA1","EOMES","KHDRBS2","RBM24","UNC5D","OAS1")
)

# ---------------------------------------------
# split dataset into train/test in group-wise manner. Note that we 
# include patients not in any particular group because they serve
# as negatives
TT_STATUS <- character(nrow(MB.pheno))
for (g in unique(MB.pheno$STATUS)) {
    idx <- which(MB.pheno$STATUS %in% g)
    status <- rep("TEST",length(idx))
    status[1:(floor(pctTrain * length(idx)))] <- "TRAIN"
    TT_STATUS[idx] <- sample(status, replace=FALSE) # scramble
}
MB.pheno <- cbind(MB.pheno, TT_STATUS=TT_STATUS)

# ---------------------------------------------
# build predictor for each subtype
MB.pheno_train <- subset(MB.pheno, TT_STATUS %in% "TRAIN")

# custom similarity measure
geneSim <- function(x) {
    if (nrow(x)>=1) x <- x[1,]
    nm <- colnames(x)
    x <- as.numeric(x)
    n <- length(x)
    rngX  <- max(x)-min(x)
    
    out <- matrix(NA,nrow=n,ncol=n);
    # weight between i and j is
    # wt(i,j) = 1 - (abs(g[i]-g[j])/(max(g)-min(g)))
    # where g is the eMB.xpression vector for each gene
    for (j in 1:n) out[,j] <- 1-(abs((x-x[j])/rngX))
    rownames(out) <- nm; colnames(out)<- nm
    out
}

# directories with group-specific predictors
predRes <- list()
for (g in names(groupSig)){
    pDir <- sprintf("%s/%s",outDir,g)
    if (file.exists(pDir)) unlink(pDir)
    dir.create(pDir)
    
    # we want each gene to have its own PSN
    sigNets <- list()
    for (g2 in groupSig[[g]]) sigNets[[g2]] <- g2
    
    # create patient networks using train & test samples
    # networks are limited to signature genes for this subtype
    idx     <- which(MB.xpr_names %in% groupSig[[g]])
    cat(sprintf("Subtype : %s { %s } => %i measures\n",
                g, paste(groupSig[[g]],collapse=","), length(idx)))
    netDir  <- sprintf("%s/networks",pDir)
    if (!file.exists(netDir)) unlink(netDir)
    netList <- makePSN_NamedMatrix(MB.xpr[idx,], MB.xpr_names[idx], 
                                   sigNets,netDir,
                                   simMetric="custom",customFunc=geneSim,
                                   verbose=TRUE)
    
    # create a GeneMANIA database out of these networks
    dbDir <- GM_createDB(netDir, MB.pheno$ID, pDir)
    
    # run a query using training samples for this subtype.
    # get ranking for all patients in the database
    trainSamps <- MB.pheno$ID[which(MB.pheno$TT_STATUS %in% "TRAIN" 
								 & MB.pheno$STATUS %in% g)]
    qFile      <- sprintf("%s/query.txt", pDir)
    GM_writeQueryFile(trainSamps, "all", nrow(MB.pheno),qFile)
    resFile <- runGeneMANIA(dbDir$dbDir, qFile, pDir)
    
    # compute ROC curve for each predictor
    predRes[[g]] <- GM_getQueryROC(sprintf("%s.PRANK",resFile),MB.pheno, g)
    
}

par(mfrow=c(2,2))
tmp <- sapply(names(predRes), function(nm){
    x <- predRes[[nm]]
    plot(x$roc,
         main=sprintf("%s: N=%i (AUC=%1.3f)",
                      nm,length(x$roc@x.values[[1]]),x$auc),
         cex.main=0.8)
    })

save(predRes,file=sprintf("%s/predictions.Rdata",outDir))

# ---------------------------------------------
# finally, predict the class of each test sample
predClass 	<- GM_OneVAll_getClass(predRes)
testSamps 	<- merge(x=MB.pheno,y=predClass,by="ID")

### compute class match accuracy
cat("\n-------------------------------------\n")
rightClass <- testSamps$STATUS == testSamps$PRED_CLASS
numCor <- sum(rightClass); ln <- nrow(testSamps)
cat(sprintf("Overall classifier accuracy = %i of %i  (%i%%)",
            numCor, ln, round((numCor/ln)*100)))
# class-specific accuracy
cat("\nClass-specific accuracy: \n")
for (g in names(groupSig)){
    idx <- which(testSamps$STATUS %in% g)
    numCor <- sum(rightClass[idx]); ln <- length(idx)
    cat(sprintf("\t%s = %i of %i  (%i%%)\n",
            g, numCor, ln, round((numCor/ln)*100)))
}
