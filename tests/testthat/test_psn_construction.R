

###test_that("invalid sims is flagged",{
###
###    expect_equal(TRUE, checkSimValid(list(a="pearsonCorr")))
###    expect_identical(TRUE, checkSimValid(list(a="pearsonCorr",b=function(x) 2+4)))
###    expect_identical(TRUE, checkSimValid(list(a="normDiff")))
###    expect_error(checkSimValid(list(a="normDifff")))
###    expect_error(checkSimValid(list(a=2)))
###})

#test_that("psns constructed using provided sims",{
    rm(list=ls())
    library(curatedTCGAData)
    library(netDx)
    brca <- suppressMessages(
    curatedTCGAData("BRCA",
               c("mRNAArray","RPPA*","Methylation_methyl27*"),
	dry.run=FALSE,version="1.1.38"))
    
    source("prepare_data.R")
    brca <- prepareData(brca,setBinary=TRUE)

    groupList <- list()

    # genes in mRNA data are grouped by pathways
    pathFile <- "pathway_ex3.gmt"
    pathList <- readPathways(pathFile)
    groupList[["BRCA_mRNAArray-20160128"]] <- pathList
    # clinical data is not grouped; each variable is its own feature
    groupList[["clinical"]] <- list(
          age="patient.age_at_initial_pathologic_diagnosis",
    	   stage="STAGE"
    )
    # for methylation generate one feature containing all probes
    # same for proteomics data
    for (k in 2:3){
    tmp <- list(rownames(experiments(brca)[[k]]));
    names(tmp) <- names(brca)[k]
    groupList[[names(brca)[k]]] <- tmp
    }

dataList <- dataList2List(brca,groupList)
pheno <- data.frame(INTERNAL_ID=1:nrow(colData(brca)),
    ID=colData(brca)$ID)

netDir <- paste(tempdir(),"nets",sep=getFileSep())
# ----------
message("test pearson and built-in function")
if (file.exists(netDir)) unlink(netDir,recursive=TRUE)
dir.create(netDir,recursive=TRUE)
sims <- list(a="pearsonCorr",b="normDiff",
        c="pearsonCorr",d="pearsonCorr")
names(sims) <- names(groupList)
x <- createNetFuncFromSimList(dataList$assays,groupList,
    netDir,sims)
prof <- dir(netDir,".profile$")
cont <- dir(netDir,"_cont.txt")
testthat::expect_equal(7,length(x))
pathways <- c(names(groupList[[1]]),names(groupList[3:4]))
testthat::expect_identical(sort(pathways), sort(sub(".profile$","",prof)))
print(all.equal(sort(names(groupList[[2]])),sort(sub("_cont.txt","",cont))))

# ----------
message("just pearson")
gp1 <- groupList[1]; d1 <- dataList$assays[1]
sims <- list(a="pearsonCorr"); names(sims) <- names(gp1)
if (file.exists(netDir)) unlink(netDir,recursive=TRUE)
dir.create(netDir,recursive=TRUE)
x <- createNetFuncFromSimList(d1,gp1,
    netDir,sims)
prof <- dir(netDir,".profile$")
cont <- dir(netDir,"_cont.txt")
testthat::expect_equal(length(prof),3)
testthat::expect_equal(length(cont),0)

message("just built-in")
gp1 <- groupList[2]; d1 <- dataList$assays["clinical"]
sims <- list(a="normDiff"); names(sims) <- names(gp1)
if (file.exists(netDir)) unlink(netDir,recursive=TRUE)
dir.create(netDir,recursive=TRUE)
x <- createNetFuncFromSimList(d1,gp1,
    netDir,sims)
prof <- dir(netDir,".profile$")
cont <- dir(netDir,"_cont.txt")
testthat::expect_equal(length(prof),0)
testthat::expect_equal(length(cont),2)

message("just custom")
gp1 <- groupList[2]; d1 <- dataList$assays["clinical"]
sims <- list(a=normDiff); names(sims) <- names(gp1)
if (file.exists(netDir)) unlink(netDir,recursive=TRUE)
dir.create(netDir,recursive=TRUE)
x <- createNetFuncFromSimList(d1,gp1,
    netDir,sims)
prof <- dir(netDir,".profile$")
cont <- dir(netDir,"_cont.txt")
testthat::expect_equal(length(prof),0)
testthat::expect_equal(length(cont),2)

message("test full net creation function")
sims <- list(a="pearsonCorr",b="normDiff",
        c="pearsonCorr",d="pearsonCorr")
names(sims) <- names(groupList)
if (file.exists(netDir)) unlink(netDir,recursive=TRUE)
dir.create(netDir,recursive=TRUE)
netList<- createPSN_MultiData(
    dataList=dataList$assays,
    groupList=groupList,
    pheno=pheno, netDir=netDir,
    makeNetFunc=NULL, sims=sims,
    numCores=1,verbose=TRUE)
expect_equal(7, length(netList))
testthat::expect_equal(2, length(dir(sprintf("%s/INTERACTIONS",netDir))))


# test full predictor
#set.seed(42) # make results reproducible
##out <- 
##   buildPredictor(dataList=brca,groupList=groupList,
##      sims=sims,
##      outDir=outDir, ## netDx requires absolute path
##      numSplits=2L, featScoreMax=2L, featSelCutoff=1L,
##	  numCores=nco)
##    
#})