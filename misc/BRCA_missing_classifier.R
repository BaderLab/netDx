# runs netdx classifier with brca {xpr,dnam,mirna}
# compares singletons to 2 datatypes to all three.
# 4 way classifier
rm(list=ls())
require(netDx)

dt <- format(Sys.Date(),"%y%m%d")
dataFile <- "/Users/shraddhapai/Documents/Research/BaderLab/TCGA_breastCancer/input/BRCA_3data_input_170117.Rdata"

outDir <- sprintf("/Users/shraddhapai/Documents/Research/BaderLab/TCGA_breastCancer/output/missing_%s", dt)

numCores <- 2L
set.seed(1203);

if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
dir.create(outDir)

logFile <- sprintf("%s/log.txt", outDir)
sink(logFile,split=TRUE)
tryCatch({

# load data
cat("loading\n")
load(dataFile)
subtypes <- unique(pheno$PAM50)
subtypes <- setdiff(subtypes,"Normal")
print(subtypes)
print(table(pheno$PAM50))
colnames(pheno)[4] <- "STATUS"


rm(dnam,mirna)

#alldat <- rbind(rbind(xpr,dnam),mirna)
# elements to group into one network per datatype
netSets <- list(xpr=rownames(xpr))

# split into train/test for all 4 classes.
pheno$TT_STATUS <- splitTestTrain(pheno)

# run prediction using each data combination in turn
N <- nrow(xpr)*ncol(xpr)
for (pctMiss in c(10,50,70,85,90,95,99)) {
	t0 <- Sys.time()
	pDir <- sprintf("%s/miss%s",outDir, pctMiss)
	dir.create(pDir)


	tmp_xpr <- xpr
	
	# insert missing values
	idx <- sample(N,floor((pctMiss/100)*N),FALSE)
	#i <- idx %% nrow(tmp_xpr); 
	#if (any(i==0)) i[which(i==0)] <- nrow(tmp_xpr)
	#j <- ceiling(idx/nrow(tmp_xpr))
	tmp_xpr <- as.matrix(tmp_xpr)
	tmp_xpr[idx] <- NA

	cat(sprintf("%% missing = %i ; %i entries set to NA\n",
			pctMiss, sum(is.na(tmp_xpr))))

	profDir <- sprintf("%s/profiles",pDir)
	dir.create(profDir)

	# create similarity networks for all datatypes
	tmp <- makePSN_NamedMatrix(tmp_xpr, rownames(xpr),netSets,
						   profDir,verbose=FALSE, numCores=numCores,
						   writeProfiles=TRUE)

	# create a single GM database with training and test samples
	dbDir 	<- GM_createDB(profDir,pheno$ID,pDir,numCores=numCores)

	outRes <- list()
	for (g in subtypes) {
		qSamps <- pheno$ID[which(pheno$STATUS %in% g & 
								 pheno$TT_STATUS%in% "TRAIN")]
		qFile <- sprintf("%s/%s_testQuery",pDir,g)

		# use only selected nets
		GM_writeQueryFile(qSamps,"all",nrow(pheno),qFile)

		resFile <- runGeneMANIA(dbDir$dbDir,qFile,resDir=pDir)
		outRes[[g]] <- GM_getQueryROC(sprintf("%s.PRANK",resFile),pheno,g)
	}
	outClass <- GM_OneVAll_getClass(outRes)
	both <- merge(x=pheno,y=outClass,by="ID")
	
	save(outRes,file=sprintf("%s/outRes.Rdata",pDir))
	write.table(both,file=sprintf("%s/predictionResults.txt",pDir),
				sep="\t",col=TRUE,row=FALSE,quote=FALSE)

	cat(sprintf("%s complete\n", pctMiss))
	print(Sys.time()-t0)
}

},error=function(ex){
	print(ex)
},finally={
	cat("Closing log\n")
	sink(NULL)
})
