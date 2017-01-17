# runs netdx classifier with brca {xpr,dnam,mirna}
# compares singletons to 2 datatypes to all three.
# 4 way classifier

dataFile <- "/Users/shraddhapai/Documents/Research/BaderLab/TCGA_breastCancer/input/BRCA_3data_input_170117.Rdata"
outDir <- "/Users/shraddhapai/Documents/Research/BaderLab/TCGA_breastCancer/output/data3_170117"
numCores <- 2L

# load data
cat("loading\n")
load(dataFile)
subtypes <- unique(pheno$PAM50)
subtypes <- setdiff(subtypes,"Normal")
print(subtypes)
print(table(pheno$PAM50))
colnames(pheno)[4] <- "STATUS"

profDir <- sprintf("%s/profiles",outDir)
dir.create(profDir)

alldat <- rbind(rbind(xpr,dnam),mirna)

# elements to group into one network per datatype
netSets <- list(
	xpr=rownames(xpr),
	dnam=rownames(dnam),
	mirna=rownames(mirna))
browser()

# create similarity networks for all datatypes.
require(netDx)
tmp <- makePSN_NamedMatrix(alldat, rownames(alldat),netSets,
						   profDir,verbose=FALSE, numCores=numCores,
						   writeProfiles=TRUE)
# create a single GM database with training and test samples
dbDir 	<- GM_createDB(profDir,pheno$ID,outDir,numCores=numCores)
# split into train/test for all 4 classes.
pheno$TT_STATUS <- splitTestTrain(pheno)

# for each combination of datatypes (start with just rna)
outRes <- list()
for (g in subtypes) {
	qSamps <- pheno$ID[which(pheno$STATUS %in% g & 
							 pheno$TT_STATUS%in% "TRAIN")]
	qFile <- sprintf("%s/%s_testQuery",outDir,g)
	GM_writeQueryFile(qSamps,"all",nrow(pheno),qFile)
	resFile <- runGeneMANIA(dbDir$dbDir,qFile,resDir=outDir)
	system(sprintf("unlink %s",resFile))
	outRes[[g]] <- GM_getQueryROC(sprintf("%s.PRANK",resFile),pheno,g)
}
outClass <- GM_OneVAll_getClass(outRes)
both <- merge(x=pheno,y=outClass,by="ID")
both <- both[-which(both$STATUS %in% "Normal"),]
acc <- sum(both$STATUS == both$PRED_CLASS)/nrow(both)
print(table(both[,c("STATUS","PRED_CLASS")]))

save(outRes,file=sprintf("%s/outRes.Rdata",outDir))
write.table(both,file=sprintf("%s/predictionResults.txt",outDir),sep="\t",
			col=TRUE,row=FALSE,quote=FALSE)


# for each class k, run gm with query of class k.
# end outer for


