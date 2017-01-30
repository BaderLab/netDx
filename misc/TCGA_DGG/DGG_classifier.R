# runs netdx classifier with diffuse grade glioma
# compares singletons to 2 datatypes 
# 3 way classifier
rm(list=ls())

dt <- format(Sys.Date(),"%y%m%d")
dataFile <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_DGG/input/DGG_input_170118.Rdata"

outDir <- sprintf("/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_DGG/output/integrate_%s", dt)

numCores <- 2L

if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
dir.create(outDir)

logFile <- sprintf("%s/log.txt", outDir)
sink(logFile,split=TRUE)
tryCatch({

# load data
cat("loading\n")
load(dataFile)
subtypes <- unique(pheno$STATUS)
print(subtypes)
print(table(pheno$STATUS))

profDir <- sprintf("%s/profiles",outDir)
dir.create(profDir)

alldat <- rbind(xpr,dnam)
# elements to group into one network per datatype
netSets <- list(
	xpr=rownames(xpr),
	dnam=rownames(dnam)
	)

# create similarity networks for all datatypes.
require(netDx)
tmp <- makePSN_NamedMatrix(alldat, rownames(alldat),netSets,
						   profDir,verbose=FALSE, numCores=numCores,
						   writeProfiles=TRUE)

# create a single GM database with training and test samples
dbDir 	<- GM_createDB(profDir,pheno$ID,outDir,numCores=numCores)

# split into train/test for all 4 classes.
pheno$TT_STATUS <- splitTestTrain(pheno)

# datatype combinations to try.
combList <- list(xpr="xpr",dnam="dnam",all=c("xpr","dnam"))

# run prediction using each data combination in turn
for (cur in names(combList)) {
	t0 <- Sys.time()
	cat(sprintf("%s\n",cur)) 
	pDir <- sprintf("%s/%s",outDir, cur)
	dir.create(pDir)

	outRes <- list()
	for (g in subtypes) {
		qSamps <- pheno$ID[which(pheno$STATUS %in% g & 
								 pheno$TT_STATUS%in% "TRAIN")]
		qFile <- sprintf("%s/%s_testQuery",pDir,g)

		# use only selected nets
		if (cur == "all") nets <- "all" 
		else nets <- paste(combList[[cur]],".profile",sep="")
		GM_writeQueryFile(qSamps,nets,nrow(pheno),qFile)

		resFile <- runGeneMANIA(dbDir$dbDir,qFile,resDir=pDir)
		outRes[[g]] <- GM_getQueryROC(sprintf("%s.PRANK",resFile),pheno,g)
	}
	outClass <- GM_OneVAll_getClass(outRes)
	both <- merge(x=pheno,y=outClass,by="ID")

	#both <- both[-which(both$STATUS %in% "Normal"),]
	acc <- sum(both$STATUS == both$PRED_CLASS)/nrow(both)
	print(table(both[,c("STATUS","PRED_CLASS")]))
	
	save(outRes,file=sprintf("%s/outRes.Rdata",pDir))
	write.table(both,file=sprintf("%s/predictionResults.txt",pDir),
				sep="\t",col=TRUE,row=FALSE,quote=FALSE)

	cat(sprintf("%s complete\n", cur))
	print(Sys.time()-t0)
}

},error=function(ex){
	print(ex)
},finally={
	cat("Closing log\n")
	sink(NULL)
})
