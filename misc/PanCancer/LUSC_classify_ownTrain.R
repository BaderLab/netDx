#' classify by survival for TCGA SKCM
#' uses own train/test split

inDir <- "/mnt/data2/BaderLab/PanCancer_LUSC/input"
outRoot <-"/mnt/data2/BaderLab/PanCancer_LUSC/output"

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/ownTrain_%s",outRoot,dt)
# ----------------------------------------------------------------
# helper functions

# normalized difference 
# x is vector of values, one per patient (e.g. ages)
normDiff <- function(x) {
    #if (nrow(x)>=1) x <- x[1,]
    nm <- colnames(x)
    x <- as.numeric(x)
    n <- length(x)
    rngX  <- max(x,na.rm=T)-min(x,na.rm=T)
    
    out <- matrix(NA,nrow=n,ncol=n);
    # weight between i and j is
    # wt(i,j) = 1 - (abs(x[i]-x[j])/(max(x)-min(x)))
    for (j in 1:n) out[,j] <- 1-(abs((x-x[j])/rngX))
    rownames(out) <- nm; colnames(out)<- nm
    out
}
# -----------------------------------------------------------
# process input
inFiles <- list(
	clinical=sprintf("%s/LUSC_clinical_core.txt",inDir),
	survival=sprintf("%s/LUSC_binary_survival.txt",inDir),
	rna=sprintf("%s/LUSC_mRNA_core.txt",inDir),
	mirna=sprintf("%s/LUSC_miRNA_core.txt",inDir),
	rppa=sprintf("%s/LUSC_RPPA_core.txt",inDir),
	train=sprintf("%s/LUSC_train_sample_list.txt",inDir),
	test=sprintf("%s/LUSC_test_sample_list.txt",inDir)
)

pheno <- read.delim(inFiles$clinical,sep="\t",h=T,as.is=T)
colnames(pheno)[1] <- "ID"

surv <- read.delim(inFiles$survival,sep="\t",h=T,as.is=T)
colnames(surv)[1:2] <- c("ID","STATUS_INT")
survStr <- rep(NA,nrow(surv))
survStr[surv$STATUS_INT<1] <- "SURVIVENO"
survStr[surv$STATUS_INT>0] <- "SURVIVEYES"
surv$STATUS <- survStr
pheno <- merge(x=pheno,y=surv,by="ID")


dats <- list()

# rnaseq 
cat("\t* RNA\n")
t0 <- Sys.time()
rna <- read.delim(inFiles$rna,sep="\t",h=T,as.is=T)
print(Sys.time()-t0)
rna <- t(rna)
colnames(rna) <- rna[1,]; rna <- rna[-1,]; 
rna <- rna[-nrow(rna),]
class(rna) <- "numeric"
dats$rna <- rna
rm(rna)

# mirna
cat("\t* miR\n")
t0 <- Sys.time()
mir <- read.delim(inFiles$mirna,sep="\t",h=T,as.is=T)
print(Sys.time()-t0)
mir <- t(mir)
colnames(mir) <- mir[1,]; mir <- mir[-1,]; 
mir <- mir[-nrow(mir),]
class(mir) <- "numeric"
dats$mir <- mir; rm(mir)

# Proteomics
cat("\t* RPPA\n")
rppa <- read.delim(inFiles$rppa,sep="\t",h=T,as.is=T)
rppa <- t(rppa)
colnames(rppa) <- rppa[1,]; rppa <- rppa[-1,]; 
rppa <- rppa[-nrow(rppa),]
class(rppa) <- "numeric"
dats$rppa <- rppa; rm(rppa)

# include only data for patients in classifier
dats <- lapply(dats, function(x) { x[,which(colnames(x)%in%pheno$ID)]})
dats <- lapply(dats, function(x) { 
	midx <- match(pheno$ID,colnames(x))
	x <- x[,midx]
	x
})

# clinical
cat("\t* Clinical\n")
clin <- pheno
# ######
# This section copied from main.R of syn1895966 and adapted to the current
# code
clin$stage <- as.vector(clin$stage)
clin$stage[clin$stage=="Stage IA"| clin$stage=="Stage IB"] <- "I"
clin$stage[clin$stage=="Stage IIA"| clin$stage=="Stage IIB"| clin$stage=="Stage II"] <- "II"
clin$stage[clin$stage=="Stage IIIA"| clin$stage=="Stage IIIB"] <- "III"
clin$stage <- as.factor(clin$stage)
clin <- clin[, -which(colnames(clin)=="gender")]
# ######
rownames(clin) <- clin[,1]; 
clin <- t(clin[,c("age","stage")])
clin[1,] <- as.integer(clin[1,])
clin[2,] <- as.integer(as.factor(clin[2,]))
class(clin) <- "numeric"
dats$clinical <- clin; rm(clin)

alldat <- do.call("rbind",dats)

cat("Sanity check, data type dims\n")
tmp <- sapply(names(dats), function(x) {
	cat(sprintf("%s\n", x)); print(dim(dats[[x]]));
})

#train <- read.delim(inFiles$train,sep="\t",h=F,as.is=T)
#test <-read.delim(inFiles$test,sep="\t",h=F,as.is=T) 

# ----------------------------------------------------------
# build classifier
numCores <- 8L
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
dir.create(outDir)

# input nets for each category
netSets <- lapply(dats, function(x) rownames(x))

# create similarity networks for all datatypes.
require(netDx)
netDir <- sprintf("%s/networks",outDir)
dir.create(netDir)
cat("about to make nets\n")
tmp <- makePSN_NamedMatrix(alldat, rownames(alldat),netSets,
	   netDir,verbose=FALSE, numCores=numCores,
	   writeProfiles=TRUE)
cat("about to make db\n")
# create a single GM database with training and test samples
dbDir 	<- GM_createDB(netDir,pheno$ID,outDir,numCores=numCores)

# datatype combinations to try.
combList <- list(
	clinical="clinical.profile",
	clinicalArna=c("clinical.profile","rna.profile"),
	clinicalAmir=c("clinical.profile","mir.profile"),
	clinicalArppa=c("clinical.profile","rppa.profile"),
	all="all")

subtypes <- unique(pheno$STATUS)


# ---------------------------------------------------------------------
# run test for different train/test splits
sink(sprintf("%s/log.txt",outDir),split=TRUE)
tryCatch({

for (runNum in 1:100) {
	curd <- sprintf("%s/run%i",outDir,runNum)
	dir.create(curd)
	cat(sprintf("Run %i\n", runNum))

	#TT_STATUS <- rep(NA,nrow(pheno))
	#TT_STATUS[which(pheno$ID %in% train[,runNum])] <- "TRAIN"
	#TT_STATUS[which(pheno$ID %in% test[,runNum])] <- "TEST"
	TT_STATUS <- splitTestTrain(pheno$ID,pctTrain=0.8,setseed=runNum*10);
	pheno$TT_STATUS <- TT_STATUS
	print(table(pheno[c("STATUS","TT_STATUS")],useNA="always"))
	is_na <- which(is.na(pheno$TT_STATUS))
	if (any(is_na)) pheno$TT_STATUS[is_na] <- "TRAIN"

	x <- sum(pheno$STATUS %in% "SURVIVENO" & pheno$TT_STATUS %in% "TRAIN")
	x <- x/sum(pheno$STATUS%in% "SURVIVENO")
	x <- round(x*100)
	y <- sum(pheno$STATUS %in% "SURVIVEYES" & pheno$TT_STATUS %in% "TRAIN")
	y <- y/sum(pheno$STATUS%in% "SURVIVEYES")
	y <- round(y*100)
	cat(sprintf("0: %i%% train ; 1: %i%% train\n",x,y) )


	# run prediction using each data combination in turn
	for (cur in names(combList)) {
		t0 <- Sys.time()
		cat(sprintf("%s\n",cur)) 
		pDir <- sprintf("%s/%s",curd, cur)
		dir.create(pDir)
	
		outRes <- list()
		for (g in subtypes) {
			qSamps <- pheno$ID[which(pheno$STATUS %in% g & 
				 pheno$TT_STATUS%in% "TRAIN")]
			qFile <- sprintf("%s/%s_testQuery",pDir,g)
	
			# use only selected nets
			GM_writeQueryFile(qSamps,combList[[cur]],nrow(pheno),qFile)
	
			resFile <- runGeneMANIA(dbDir$dbDir,qFile,resDir=pDir)
			outRes[[g]] <- GM_getQueryROC(
				sprintf("%s.PRANK",resFile),pheno,g)
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
}

},error=function(ex){
	print(ex)
},finally={
	sink(NULL)
})


