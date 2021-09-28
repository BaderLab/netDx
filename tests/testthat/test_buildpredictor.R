# test buildPredictor.R components

# ------------------------------------
# setup 

# load libraries
options(stringsAsFactors = FALSE)

# ---------------------------------------------
# tests start

test_that("split test/train works", {
	set.seed(42)
	pct <- 0.8
	pheno <- data.frame(ID=sprintf("PAT%i",1:20),
		STATUS=rep(c("case","control"),each=10))

	TT_STATUS <- splitTestTrain(pheno,pctT=0.8,verbose=FALSE)
	expect_equal(length(TT_STATUS),nrow(pheno))
	# pct matches that requested
	expect_equal(sum(TT_STATUS == "TRAIN"),16)
	# split is even among classes
	expect_equal(sum(TT_STATUS=="TRAIN" & pheno$STATUS == "case"),pct*10)
	expect_equal(sum(TT_STATUS=="TRAIN" & pheno$STATUS == "control"),pct*10)
	expect_equal(sum(TT_STATUS=="TEST" & pheno$STATUS == "case"),(1-pct)*10)

})

test_that("feature construction and compilation works",{
	# 20 patients, 10 case, 10 control
	pheno <- data.frame(ID=sprintf("PAT%i",1:20),
		STATUS=rep(c("case","control"),each=10))
	# 100 dummy genes
	rna <- matrix(rnorm(100*20),nrow=100); 
	colnames(rna) <- pheno$ID
	rownames(rna) <- sprintf("gene%i",1:100)	
	# 2 dummy clin variables
	clin <- t(data.frame(AGE=runif(20,10,50)))
	colnames(clin) <- pheno$ID

	# netDx files
	dataList <- list(rna=rna,clinical=clin)
	groupList <-list(rna=list(group1=sprintf("gene%i",seq_len(10)),
														group2=sprintf("gene%i",seq(2,100,2))),
								   clinical=list(age=c("AGE")))

 makeNets <- function(dataList, groupList, netDir,...) {
    netList <- c()
    # make RNA nets: group by pathway
    netList <- makePSN_NamedMatrix(dataList[["rna"]],
                rownames(dataList[["rna"]]),
                groupList[["rna"]],
                netDir,verbose=FALSE,
                writeProfiles=TRUE,...)
    netList <- unlist(netList)
    message(sprintf("Made %i RNA pathway nets\n", length(netList)))

    # make clinical nets,one net for each variable
    netList2 <- makePSN_NamedMatrix(dataList$clinical,
        rownames(dataList$clinical),
        groupList[["clinical"]],netDir,
        simMetric="custom",customFunc=normDiff, # custom function
        writeProfiles=FALSE,
        sparsify=TRUE,verbose=TRUE,...)
    netList2 <- unlist(netList2)
    netList <- c(netList,netList2)
    return(netList)
}
	# directory contains GENES.TXT, NETWORKS.TXT INTERACTIONS folder
	 outDir <- tempdir()
	netDir <- sprintf("%s/tmp",outDir)
	if (file.exists(netDir)) unlink(netDir,recursive=TRUE)
	dir.create(netDir)

	pheno_id <- setupFeatureDB(pheno,netDir)
	x <-createPSN_MultiData(dataList=dataList,groupList=groupList,
			pheno=pheno_id,
			netDir=netDir,makeNetFunc=makeNets,numCores=1,
			verbose=FALSE)
	
	# number of nets equals those submitted for creation
	y <- c("group1.profile","group2.profile","age_cont.txt")
	expect_equal(sort(x),sort(y))

	# created in dir
	#expect_equal(dir(sprintf("%s/INTERACTIONS",netDir),pattern="cont.txt"),
	#	"age_cont.txt")
	#expect_equal(length(dir(sprintf("%s/profiles",netDir),pattern="profile$")),,length(grep("profile$",y)))

	# custom similarity function was used as provided
	tmp <- read.delim(sprintf("%s/INTERACTIONS/1.3.txt",netDir),sep="\t",
			header=FALSE,as.is=TRUE)
	x1 <- tmp[1,1]; x2 <- tmp[1,2]
	z <- normDiff(clin)
	expect_equal(round(z[x1,x2],3),round(tmp[1,3],3))

	# compiling features 
	dbDir <- compileFeatures(netDir,outDir, numCores=1,verbose=TRUE)
	oDir <- dbDir$netDir
	expect_equal(length(dir(sprintf("%s/INTERACTIONS",oDir), ".txt")),3) # three networks
	expect_match(dir(oDir, "GENES.txt"),"GENES.txt")
	expect_match(dir(oDir, "NETWORKS.txt"),"NETWORKS.txt")
	expect_match(dir(dbDir$dbDir,"lucene.index"),"lucene.index") # copied ok
})

test_that("holding out and validation set prediction works as it should", {

	data(xpr,pheno)
	pheno$AGE <- runif(nrow(pheno),min=10,max=45)
	pathList <- list(
		a=sample(rownames(xpr),25,F),
		b=sample(rownames(xpr),100,F),
		c=sample(rownames(xpr),30,F)
	)
	suppressMessages(require(MultiAssayExperiment))
	objlist <- list(a=SummarizedExperiment(xpr))
	mae <- MultiAssayExperiment(objlist,pheno)
	
	# test validation set features
	dset <- subsampleValidationData(mae,pctValidation=0.2)
	testthat::expect_equal(nrow(colData(dset[[1]])) + nrow(colData(dset[[2]])),
				nrow(colData(mae)) )
	tst <- colData(dset[[2]])$STATUS
	lbl <- unique(tst)
	testthat::expect_equal(sum(tst == lbl[1]), sum(tst==lbl[2]))
	testthat::expect_equal(0.2, round(length(tst)/nrow(colData(mae)),digits=1) )

	outDir <- tempdir()
	if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
	dir.create(outDir)
	featSel <- list(LumA=c("a.profile","age_cont.txt"),notLumA=c("a.profile","b.profile"))
	predModel <- suppressMessages(
  		predict(trainMAE=dset[[1]], testMAE=dset[[2]], 
    		groupList=list(a=pathList,clinical=list(age="AGE")), 
    		selectedFeatures=featSel,
    		sims=list(a="pearsonCorr",clinical="normDiff"),
    		outDir=outDir, verbose = FALSE)
	)

	for (g in c("LumA","notLumA")){
		luma <- read.delim(sprintf("%s/%s/%s_query",outDir,g,g),sep=":")
		pat <- unlist(strsplit(trimws(luma[1,1]),"\t")) 
		# check that query patients are training samples for the corresponding class
		testthat::expect_equal(sort(pat),sort(colData(dset[[1]])$ID[which(colData(dset[[1]])$STATUS==g)]))
		# check that only selected features are used
		nets <- read.delim(sprintf("%s/%s/networks/NETWORKS.txt",outDir,g),sep="\t",h=F)
		testthat::expect_identical(featSel[[g]],nets[,2])
		# test that only and all holdout set patients are being classified
		dat <- read.delim(sprintf("%s/%s/%s_query-results.report.txt.PRANK",outDir,g,g),h=T,skip=1,sep="\t")
		sb <- dat[which(!is.na(dat[,2])),]
		testthat::expect_equal(nrow(colData(dset[[2]])), length(intersect(sb[,1],colData(dset[[2]])$ID)))		
	}

})

